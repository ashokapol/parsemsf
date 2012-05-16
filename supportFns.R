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

protIDsplitter <- function(cur.str) {
    #browser()
    pos <- regexpr("gi[0-9]+",cur.str)
    if (pos[1] == 1)
      cur.str <- gsub("gi","gi\\|", cur.str)
    pos <- regexpr("gi\\|[0-9]+",cur.str)
    #if (pos[1] == -1)
    #  pos <- regexpr("gi[0-9]+",cur.str)
    if (pos[1] == -1)
      pos <- regexpr("IPI\\:IPI[0-9]+",cur.str)
    if (pos[1] == -1)
      pos <- regexpr("[a-z]+\\|[a-zA-Z0-9]+",cur.str) #Uniprot
    if (pos[1] != -1){
      init <- pos[1]
      last <- attr(pos, "match.length")
      cur.ID <- substr(cur.str, init, last)
      if (last+2 < nchar(cur.str))
        cur.Ann <- substr(cur.str, last+2, nchar(cur.str))
      else
        cur.Ann <- "No info found"
    }
    else {
      x <- unlist(strsplit(cur.str,"\\|"))
      cur.ID <- x[1]
      if (length(x) > 1)
        cur.Ann <- x[2]
      else
        cur.Ann <- "No info found"
    }
    return(c(cur.ID, cur.Ann));
}


##########################################################################################
updateSeq <- function(x)
{
    #PeptideID, Sequence, Position, ModificationName, TerminalModificationName, ConfidenceLevel
    SeqOrg <- as.character(x[2])
    Seq <- SeqOrg

    posn <- as.numeric(unlist(strsplit(as.character(x[3]),",")))
    idx <- order(posn)

    modn <- unlist(strsplit(as.character(x[4]),","))
    tmod <- as.character(x[6])

    confLs <- c("Low","Medium","High")
    confL <- confLs[as.numeric(x[5])]

    posn <- posn[idx]
    modn <- modn[idx]

    Np <- length(posn)
    Nm <- length(modn)

    mod <- character(0)
    if (!is.na(tmod))
    {
        mod <- paste("N-Term(", tmod, ")", sep="")
    }

    if (Np == Nm)
    {
        for (i in 1:Np)
        {
            #browser()
            N <- as.numeric(posn[i]) + 1 #original index starts from 0, make it 1
            aa <- substr(Seq,N,N)
            modPre <- paste(aa, as.character(N),sep="")
            mod <- c(mod, paste(modPre, "(", modn[i], ")", sep=""))
            substr(Seq,N,N) <- tolower(aa)
        }
    }
    return(c(PeptideID=as.numeric(x[1]), Sequence=SeqOrg, ModSequence=Seq,
        Modifications=paste(mod, collapse=";"), Confidence=confL))
}