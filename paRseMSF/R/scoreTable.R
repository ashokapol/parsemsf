scoreTable <-
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
    xcID <- as.numeric(dbGetQuery(con,"SELECT ScoreID FROM ProcessingNodeScores WHERE ScoreName='XCorr';"))
    spID <- as.numeric(dbGetQuery(con,"SELECT ScoreID FROM ProcessingNodeScores WHERE ScoreName='SpScore';"))
    prID <- as.numeric(dbGetQuery(con,"SELECT ScoreID FROM ProcessingNodeScores WHERE ScoreName='ProbabilityScore';"))
    sql_scores <- paste("SELECT xc.PeptideID,
                          Pep.MatchedIonsCount ||\"/\"|| Pep.TotalIonsCount as IonsMatched,
                          xc.ScoreValue XCorr, sp.ScoreValue Sp, Prob.ScoreValue Probability
                          FROM PeptideScores AS xc, PeptideScores AS sp, PeptideScores AS Prob, Peptides AS Pep
                          WHERE xc.ScoreID=",xcID, "AND sp.ScoreID=", spID, "AND Prob.ScoreID=", prID,
                          "AND Pep.ConfidenceLevel >=", confidence,
                          "AND xc.PeptideID= sp.PeptideID AND xc.PeptideID=Prob.PeptideID
                          AND sp.PeptideID=Prob.PeptideID
                          AND Pep.PeptideID=xc.PeptideID;")
    rs_scores <- dbSendQuery(con, sql_scores)
    scores <- fetch(rs_scores, n=-1)
    dbClearResult(rs_scores);
    
    ModsTable <- pepModsTable(msfFile, minConf=minConf)
    pepScoreTable <- sqldf("SELECT p.*, sc.IonsMatched, sc.XCorr, sc.Sp, sc.Probability FROM ModsTable p,
                  scores sc WHERE p.PeptideID=sc.PeptideID GROUP BY p.PeptideID;")
    pepScoreTable <- pepScoreTable[order(as.numeric(pepScoreTable[,1])),]
    dbDisconnect(con);
    return(pepScoreTable)
}

