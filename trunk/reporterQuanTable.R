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

reporterQuanTable <- function(msfFile, minConf="High")
{
    require(RSQLite)
    require(sqldf)
    driver = dbDriver("SQLite")
    con <- dbConnect(driver, dbname=msfFile)

    confidence <- switch(minConf,
              High = 3,
              Medium = 2,
              Low = 1,
              3)

    # Proteins (ProteinID and Description)
    sql_repQuan <- paste("SELECT p.PeptideID, p.Sequence, sh.ScanEventID, riq.Mass, riq.Height
                          FROM Peptides AS p, ReporterIonQuanResults AS riq,
                          ReporterIonQuanResultsSearchSpectra AS riqss, SpectrumHeaders AS sh
                          WHERE Peptides.SpectrumID = ReporterIonQuanResultsSearchSpectra.SearchSpectrumID
                          AND Peptides.SpectrumID = SpectrumHeaders.SpectrumID
                          AND ReporterIonQuanResultsSearchSpectra.SpectrumID = ReporterIonQuanResults.SpectrumID
                          AND ScanEventID = 1 AND Peptides.ConfidenceLevel >=", confidence,";")
    rs_repQuan <- dbSendQuery(con, sql_repQuan)
    riI <- fetch(rs_repQuan, n=-1)
    dbClearResult(rs_repQuan);
    browser()
    ModsTable <- pepModsTable(msfFile, minConf=minConf)
    pepScoreTable <- sqldf("SELECT p.*, sc.IonsMatched, sc.XCorr, sc.Sp, sc.Probability FROM ModsTable p,
                  scores sc WHERE p.PeptideID=sc.PeptideID GROUP BY p.PeptideID;")

    pepScoreTable <- pepScoreTable[order(as.numeric(pepScoreTable[,1])),]

    #Cleanup
    dbDisconnect(con);
    return(pepScoreTable)
}

# Run
# M <- scoreTable(msfFile="D:/Research/R/msfRead/test_dimethyl.msf", minConf="High")