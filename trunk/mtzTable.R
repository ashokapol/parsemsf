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

mtzTable <- function(msfFile, minConf="High")
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
    sql_mtz <- paste("SELECT p.PeptideID, sh.ScanNumbers, sh.RetentionTime, sh.Mass, sh.Charge
                          FROM Peptides AS p, SpectrumHeaders AS sh
                          WHERE p.SpectrumID=sh.SpectrumID AND",
                          "p.ConfidenceLevel >=", confidence,
                          "GROUP BY p.PeptideID;")
    rs_mtz <- dbSendQuery(con, sql_mtz)
    mtzTab <- fetch(rs_mtz, n=-1)
    dbClearResult(rs_mtz);

    pepTable <- scoreTable(msfFile, minConf=minConf)
    
    mtzTable <- sqldf("SELECT p.*, mtz.ScanNumbers, mtz.RetentionTime, mtz.Mass, mtz.Charge
                          FROM pepTable p, mtzTab mtz
                          WHERE p.PeptideID=mtz.PeptideID
                          GROUP BY p.PeptideID;")

    mtzTable <- mtzTable[order(as.numeric(mtzTable[,1])),]

    #Cleanup
    dbDisconnect(con);
    return(mtzTable)
}

# Run
# M <- mtzTable(msfFile="D:/Research/R/msfRead/test_dimethyl.msf", minConf="High")