# Written by Ashoka D. Polpitiya
# for the Translational Genomics Research Institute (TGen, Phoenix, AZ)
# Copyright 2012, Translational Genomics Research Institute
# E-mail: ashoka@tgen.org
# Website: http://paRseMSF.googlecode.com
# -------------------------------------------------------------------------
#
# Licensed under the Apache License, Version 2.0; you may not use this file except
# in compliance with the License.  You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0
#
# -------------------------------------------------------------------------

getArea <- function(currData, idx)
{
  # Calculate areas based on the two most intense isotopic peak areas.
  # For single peaks NA is returned.
  outArea <- numeric(0)
  
  if (sum(idx)>0)
  {
    featureData <- currData[idx,] # Light, Medium, or Heavy
    confLs <- featureData$ConfidenceLevel
    # Select the ones with highest conf
    featureData <- featureData[confLs==max(confLs),]

    if (dim(featureData)[1]>1)
    {
      nPeps <- unique(featureData$PeptideID)
      for (k in nPeps)
      {
        Data <- featureData[featureData$PeptideID==k,]
        areas <- sort(Data$Area, decreasing=TRUE)
        if (length(areas)>1)
          outArea <- c(outArea, sum(areas[c(1,2)]))
      }
      outArea <- sort(outArea, decreasing=TRUE)
      outArea <- outArea[1]
    }
    else
      outArea <- NA
  }
  else
    outArea <- NA

  return(outArea)
}

#-----------------------------------------
calcDimethylRatios <- function(currData, Seq, minConf=3)
{
  #if (Seq=="APNHAVVTRK")
  #  browser()
  if (max(currData$ConfidenceLevel) >= minConf) # at least one peptide is "High" confident
  {
    lightIdx <- (currData$QuanChannelID==1)
    mediumIdx <- (currData$QuanChannelID==2)
    heavyIdx <- (currData$QuanChannelID==3)
    
    Light <- getArea(currData,lightIdx)
    Medium <- getArea(currData,mediumIdx)
    Heavy <- getArea(currData,heavyIdx)

    Medium.Light <- ifelse((!is.na(Medium) & !is.na(Light)), Medium/Light, NA)
    Heavy.Light <- ifelse((!is.na(Heavy) & !is.na(Light)), Heavy/Light, NA)
    #print(Seq)

    outTab <- data.frame(Sequence=Seq, MediumToLight=Medium.Light,
                HeavyToLight=Heavy.Light, LightArea=Light,
                MediumArea=Medium, HeavyArea=Heavy)
  }
  else
    outTab <- NA
    
  return(outTab)
}


#------------------------------------------------------------------------------
dimethylRatioTable <- function(msfFile, minConf="High", ratioMinConf="Low")
{
    require(RSQLite)
    require(sqldf)
    driver = dbDriver("SQLite")
    con <- dbConnect(driver, dbname=msfFile)

    if (!checkPara(msfFile, "Dimethyl"))
      stop("This does not seem to be a Dimethyl run.")

    confidence <- switch(minConf,
                    High = 3,
                    Medium = 2,
                    Low = 1,
                    3)

    rConf <- switch(ratioMinConf,
                  High = 3,
                  Medium = 2,
                  Low = 1,
                  3)

    ratioTable <- NULL

    pepTable <- mtzTable(msfFile, minConf=minConf)
    
    pepSeqs <- pepTable$Sequence
    for (k in 1:length(pepSeqs))
    {
      Seq <- pepSeqs[k]
      pepID <- pepTable$PeptideID[k]
      
      sql_rat <- paste("SELECT p.PeptideID, p.Sequence, p.ConfidenceLevel,
                      eva.IsotopePatternID, eva.QuanResultID, eva.QuanChannelID, e.Area
                      FROM Peptides p, PrecursorIonQuanResultsSearchSpectra piqrss,
                      EventAnnotations eva, Events e
                      WHERE  piqrss.SearchSpectrumID = p.SpectrumID AND
                      eva.QuanResultID = piqrss.QuanResultID AND
                      e.EventID = eva.EventID AND
                      p.Sequence=", shQuote(Seq), " AND p.ConfidenceLevel >=",
                      rConf, ";")
      rs_rat <- dbSendQuery(con, sql_rat)
      ratTab <- fetch(rs_rat, n=-1)
      ratios <- calcDimethylRatios(ratTab, Seq, confidence)
      ratioTable <- rbind(ratioTable, ratios)
    }
    dbClearResult(rs_rat);
#
#    ratioTab <- calcDimethylRatios(aucTab, confidence)
#
    #pepTable <- mtzTable(msfFile, minConf=minConf)

    aucTab2 <- sqldf("SELECT p.*, ratios.MediumToLight, ratios.HeavyToLight,
                        ratios.LightArea, ratios.MediumArea, ratios.HeavyArea
                          FROM pepTable p, ratioTable ratios
                          WHERE p.Sequence=ratios.Sequence
                          GROUP BY p.PeptideID;")

    aucTab2 <- aucTab2[order(as.numeric(aucTab2[,1])),]

    #Cleanup
    dbDisconnect(con);
    return(aucTab2)
}

# Run
msfFileN = "C:/Ashoka/Research/R/paRseMSF/SYR_PIR355A16PA1_12AUG06_H1223 Lavage_H1218 Lavage Run 2.msf"
M <- dimethylRatioTable(msfFile=msfFileN, minConf="High")
# M <- dimethylRatioTable(msfFile=File, minConf="High")
