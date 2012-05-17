pepModsTable <-
function(msfFile, minConf="High")
{
    require(RSQLite)
    driver = dbDriver("SQLite")
    con <- dbConnect(driver, dbname=msfFile)
    
    confidence <- switch(minConf,
              High = 3,
              Medium = 2,
              Low = 1, 
              3)
    sql_regMods <- paste("SELECT paam.PeptideID, p.Sequence,",
                    "group_concat(paam.Position) AS Position,",
                    "group_concat(aam.ModificationName) AS ModificationName,",
                    "p.ConfidenceLevel",
                    "FROM PeptidesAminoAcidModifications paam, AminoAcidModifications aam, Peptides p",
                    "WHERE paam.AminoAcidModificationID = aam.AminoAcidModificationID AND",
                    "paam.PeptideID = p.PeptideID AND",
                    "p.ConfidenceLevel >=", confidence,
                    "GROUP BY paam.PeptideID",
                    "ORDER BY paam.PeptideID;")                    
    rs_regMods <- dbSendQuery(con, sql_regMods)              
    
    regMods <- fetch(rs_regMods, n=-1) # fetch 3 rows from result (use -1 to fetch all rows)
    dbClearResult(rs_regMods);
    
    sql_termMods <- paste("SELECT ptm.PeptideID,",
                      "group_concat(aam.ModificationName) AS TerminalModificationName",
                      "FROM PeptidesTerminalModifications ptm, AminoAcidModifications aam, Peptides p",
                      "WHERE ptm.TerminalModificationID = aam.AminoAcidModificationID AND",
                      "ptm.PeptideID = p.PeptideID AND",
                      "p.ConfidenceLevel >=", confidence,
                      "GROUP BY ptm.PeptideID",
                      "ORDER BY ptm.PeptideID;")
    rs_termMods <- dbSendQuery(con, sql_termMods)
    
    termMods <- fetch(rs_termMods, n=-1)
    dbClearResult(rs_termMods);
    allMods <- merge(regMods, termMods, by="PeptideID", all.x=TRUE)
    ModsTable <- data.frame(t(apply(allMods,1,updateSeq)))
    
    dbDisconnect(con);
    return(ModsTable)
}

