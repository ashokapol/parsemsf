aucTable <-
function(msfFile, minConf="High")
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
    sql_auc <- paste("SELECT p.PeptideID, sum(ev.Area) AS Area
                          FROM Peptides p, Events ev, EventAreaAnnotations eva, PrecursorIonAreaSearchSpectra pis
                          WHERE ev.EventID = eva.EventID
                          AND eva.QuanResultID = pis.QuanResultID
                          AND pis.SearchSpectrumID = p.SpectrumID AND",
                          "p.ConfidenceLevel >=", confidence,
                          "GROUP BY p.PeptideID;")
    rs_auc <- dbSendQuery(con, sql_auc)
    aucTab <- fetch(rs_auc, n=-1)
    dbClearResult(rs_auc);
    pepTable <- mtzTable(msfFile, minConf=minConf)
    aucTab2 <- sqldf("SELECT p.*, auc.Area
                          FROM pepTable p, aucTab auc
                          WHERE p.PeptideID=auc.PeptideID
                          GROUP BY p.PeptideID;")
    aucTab2 <- aucTab2[order(as.numeric(aucTab2[,1])),]
    dbDisconnect(con);
    return(aucTab2)
}

