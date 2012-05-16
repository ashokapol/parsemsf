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

protAnnotations <- function(msfFile, minConf="High")
{
    require(RSQLite)
    driver = dbDriver("SQLite")
    con <- dbConnect(driver, dbname=msfFile)

    confidence <- switch(minConf,
              High = 3,
              Medium = 2,
              Low = 1,
              3)
    # Proteins (ProteinID and Description)
    sql_prots <- paste("SELECT DISTINCT pr.ProteinID AS Proteins, REPLACE(prAn.Description,'>','') AS Description",
                    "FROM Peptides p, PeptidesProteins pr, ProteinAnnotations prAn",
                    "WHERE p.PeptideID = pr.PeptideID AND",
                    "Proteins = prAn.ProteinID AND",
                    "p.ConfidenceLevel >=", confidence, ";")
    rs_prots <- dbSendQuery(con, sql_prots)
    prots <- fetch(rs_prots, n=-1)
    dbClearResult(rs_prots);
    # Split Protein Description into Protein Name and Annotation fields
    X <- matrix(unlist(lapply(prots[,2],protIDsplitter)),ncol=2,byrow=TRUE)
    # Proteins (ProteinID, ProteinName, Annotation)
    prots <- data.frame(ProteinID=prots[,1], Name=X[,1], Annotation=X[,2])

    # Peptides and Proteins
    sql_pepprots <- paste("SELECT p.PeptideID, p.Sequence, pr.ProteinID",
                      "FROM Peptides AS p LEFT OUTER JOIN PeptidesProteins AS pr",
                      "ON (p.PeptideID = pr.PeptideID)",
                      "WHERE p.ConfidenceLevel >=", confidence,
                      "ORDER BY p.PeptideID;")
    rs_pepprots <- dbSendQuery(con, sql_pepprots)
    pepsProts <- fetch(rs_pepprots, n=-1)
    dbClearResult(rs_pepprots);
    # Create a table with peptides and protein info
    pepsProts <- merge(pepsProts, prots, by="ProteinID", all.x=TRUE)

    #Cleanup
    dbDisconnect(con);
    return(list(Annotations=prots, PeptidesProteins=pepsProts))
}

# Run
# M <- protAnnotations(msfFile="D:/Research/R/msfRead/test_dimethyl.msf", minConf="High")