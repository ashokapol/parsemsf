#Files = c("aucTable.R",
#          "dimethylRatioTable.R",
#          "mtzTable.R",
#          "pepModsTable.R",
#          "protAnnotations.R",
#          "protModsTable.R",
#          "scoreTable.R",
#          "supportFns.R")
#source(Files)
rm(list=ls())
sourceDir <- function(path, trace = TRUE, ...)
{
    for (nm in list.files(path, pattern = "\\.[Rr]$")) {
       if(trace) cat(nm,":")
       source(file.path(path, nm), ...)
       if(trace) cat("\n")
    }
}
sourceDir("D:/Research/R/paRseMSF")

#M <- dimethylRatioTable(msfFile="D:/Research/R/paRseMSF/SYR_PIR355A16PA1_12AUG06_H1223 Lavage_H1218 Lavage Run 2.msf", 
#          minConf="High", ratioMinConf="Medium")
          
#M2 <- protModsTable(msfFile="C:/Ashoka/Research/R/paRseMSF/SYR_PIR355A16PA1_12AUG06_H1223 Lavage_H1218 Lavage Run 2.msf", 
#          minConf="High")
          
          
          