protIDsplitter <-
function(cur.str) {
    pos <- regexpr("gi[0-9]+",cur.str)
    if (pos[1] == 1)
      cur.str <- gsub("gi","gi\\|", cur.str)
    pos <- regexpr("gi\\|[0-9]+",cur.str)
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

