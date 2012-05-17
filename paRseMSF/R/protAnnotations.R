protAnnotations <-
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
    sql_prots <- paste("SELECT DISTINCT pr.ProteinID AS Proteins, REPLACE(prAn.Description,'>','') AS Description",
                    "FROM Peptides p, PeptidesProteins pr, ProteinAnnotations prAn",
                    "WHERE p.PeptideID = pr.PeptideID AND",
                    "Proteins = prAn.ProteinID AND",
                    "p.ConfidenceLevel >=", confidence, ";")
    rs_prots <- dbSendQuery(con, sql_prots)
    prots <- fetch(rs_prots, n=-1)
    dbClearResult(rs_prots);
    X <- matrix(unlist(lapply(prots[,2],protIDsplitter)),ncol=2,byrow=TRUE)
    prots <- data.frame(ProteinID=prots[,1], Name=X[,1], Annotation=X[,2])
    sql_pepprots <- paste("SELECT p.PeptideID, p.Sequence, pr.ProteinID",
                      "FROM Peptides AS p LEFT OUTER JOIN PeptidesProteins AS pr",
                      "ON (p.PeptideID = pr.PeptideID)",
                      "WHERE p.ConfidenceLevel >=", confidence,
                      "ORDER BY p.PeptideID;")
    rs_pepprots <- dbSendQuery(con, sql_pepprots)
    pepsProts <- fetch(rs_pepprots, n=-1)
    dbClearResult(rs_pepprots);
    pepsProts <- merge(pepsProts, prots, by="ProteinID", all.x=TRUE)
    dbDisconnect(con);
    return(list(Annotations=prots, PeptidesProteins=pepsProts))
}

