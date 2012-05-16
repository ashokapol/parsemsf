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

aucTable <- function(msfFile, minConf="High")
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
    sql_auc <- paste("SELECT p.PeptideID, sum(ev.Area) AS Area
                          FROM Peptides p, Events ev, EventAreaAnnotations eva, PrecursorIonAreaSearchSpectra pis
                          WHERE ev.EventID = eva.EventID
                          AND eva.QuanResultID = pis.QuanResultID
                          AND pis.SearchSpectrumID = p.SpectrumID AND",
                          "p.ConfidenceLevel >=", confidence,
                          "GROUP BY p.PeptideID;")
    #browser()
    rs_auc <- dbSendQuery(con, sql_auc)
    aucTab <- fetch(rs_auc, n=-1)
    dbClearResult(rs_auc);

    pepTable <- mtzTable(msfFile, minConf=minConf)

    aucTab2 <- sqldf("SELECT p.*, auc.Area
                          FROM pepTable p, aucTab auc
                          WHERE p.PeptideID=auc.PeptideID
                          GROUP BY p.PeptideID;")

    aucTab2 <- aucTab2[order(as.numeric(aucTab2[,1])),]

    #Cleanup
    dbDisconnect(con);
    return(aucTab2)
}

# Run
# M <- aucTable(msfFile="D:/Research/R/msfRead/test_dimethyl.msf", minConf="High")