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

getAreas <- function(currData)
{ 
  # Calculate areas based on the two most intense isotopic peak areas.
  # For single peaks NA is returned.
  if (dim(currData)[1]>1)
  {
    areas <- sort(currData$Area, decreasing=TRUE)
    sumArea <- sum(areas[c(1,2)])
  }
  else
    sumArea <- NA
  return(sumArea)
}

#-----------------------------------------
calcDimethylRatios <- function(tab, minConf=3)
{
  # remove any non unique entries, except the peptideID column considered.
  idx <- !duplicated(tab[,-1])
  tab <- tab[idx,]
  
  # Unique sequences
  pepSeqs <- tab$Sequence
  uSeq <- unique(pepSeqs)
  
  pID <- integer(0)
  pSequence <- character(0)
  Medium.LightR <- numeric(0)
  Heavy.LightR <- numeric(0)
  LightArea <- numeric(0)
  MediumArea <- numeric(0)
  HeavyArea <- numeric(0)
  
  for (k in 1:length(uSeq))
  {
    currSeq <- uSeq[k]
    idx <- which(pepSeqs == currSeq)
    currData <- tab[idx,]
    
    if (max(currData$ConfidenceLevel) >= minConf) # at least one peptide is "High" confident 
    {
      while (length(unique(currData$QuanResultID)) > 1) # Multiple quan results
      {
        z <- table(currData$QuanResultID)
        quanIDs <- names(which(z==max(z))) # Find one that has most isotopic peaks
        quanResult <- as.integer(quanIDs[1])
        currData_new <- currData[currData$QuanResultID==quanResult,]
        if (max(currData_new$ConfidenceLevel) < minConf) # no "High" confident ones?
          currData <- currData[-currData$QuanResultID==quanResult,]
        else
        {
          currData <- currData_new
          break; # Found "High" peptides with most istopic peaks!
        }
      }
      
      uPepID <- unique(currData$PeptideID)
      if (length(uPepID)>1)
        for (j in 1:length(uPepID))
        {
          currData_new <- currData[which(currData$PeptideID==uPepID[j]),] 
          if (max(currData_new$ConfidenceLevel) >= minConf)
          {
            currData <- currData_new
            break; # Found "High" peptides with most istopic peaks! 
          }
        }
      
      lightIdx <- (currData$QuanChannelID==1)
      mediumIdx <- (currData$QuanChannelID==2)
      heavyIdx <- (currData$QuanChannelID==3)
      
      Light <- getAreas(currData[lightIdx,])
      Medium <- getAreas(currData[mediumIdx,])
      Heavy <- getAreas(currData[heavyIdx,])
      
      if (!(is.na(Light) & is.na(Medium) & is.na(Heavy)))
      {
        Medium.Light <- ifelse((!is.na(Medium) & !is.na(Light)), Medium/Light, NA)
        Heavy.Light <- ifelse((!is.na(Heavy) & !is.na(Light)), Heavy/Light, NA)
        
        gData <- currData[lightIdx|mediumIdx|heavyIdx,]
        pID <- c(pID, gData$PeptideID[1])
        pSequence <- c(pSequence, gData$Sequence[1])
        LightArea <- c(LightArea, Light)
        MediumArea <- c(MediumArea, Medium)
        HeavyArea <- c(HeavyArea, Heavy)
        Medium.LightR <- c(Medium.LightR, Medium.Light)
        Heavy.LightR <- c(Heavy.LightR, Heavy.Light)
      }
    }
  }
  outTab <- data.frame(PeptideID=pID, Sequence=pSequence,
               MediumToLight=Medium.LightR, HeavyToLight=Heavy.LightR,
               LightArea=LightArea, MediumArea=MediumArea, HeavyArea=HeavyArea)
  return(outTab)
}


#------------------------------------------------------------------------------
dimethylRatioTable <- function(msfFile, minConf="High", ratioMinConf="Medium")
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
    # Proteins (ProteinID and Description)
    sql_auc <- paste("SELECT p.PeptideID, p.Sequence, p.ConfidenceLevel, 
                      eva.IsotopePatternID, eva.QuanResultID, eva.QuanChannelID, e.Area
                      FROM Peptides p, PrecursorIonQuanResultsSearchSpectra piqrss,
                      EventAnnotations eva, Events e
                      WHERE  piqrss.SearchSpectrumID = p.SpectrumID AND
                      eva.QuanResultID = piqrss.QuanResultID AND
                      e.EventID = eva.EventID AND",
                      "p.ConfidenceLevel >=", rConf, ";")
    
    rs_auc <- dbSendQuery(con, sql_auc)
    aucTab <- fetch(rs_auc, n=-1)
    dbClearResult(rs_auc);
    
    ratioTab <- calcDimethylRatios(aucTab, confidence)
    
    pepTable <- mtzTable(msfFile, minConf=minConf)

    aucTab2 <- sqldf("SELECT p.*, ratios.MediumToLight, ratios.HeavyToLight,
                        ratios.LightArea, ratios.MediumArea, ratios.HeavyArea
                          FROM pepTable p, ratioTab ratios
                          WHERE p.PeptideID=ratios.PeptideID;")
                          #GROUP BY p.PeptideID;")
                          
    aucTab2 <- aucTab2[order(as.numeric(aucTab2[,1])),]

    #Cleanup
    dbDisconnect(con);
    return(aucTab2)
}

# Run
# M <- dimethylRatioTable(msfFile="D:/Research/R/msfRead/test_dimethyl.msf", minConf="High")
# M <- dimethylRatioTable(msfFile="D:/Research/R/paRseMSF/SYR_PIR355A16PA1_12AUG06_H1223 Lavage_H1218 Lavage Run 2.msf", minConf="High")
