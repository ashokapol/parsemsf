updateSeq <-
function(x)
{
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

