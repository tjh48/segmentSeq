coord2gff3 <- function(cD, file, locID = "locus")
  {
    write("##gff-version 3", file)
    return(write.table(cbind(as.character(seqnames(cD@coordinates)), "rtracklayer", locID, start(cD@coordinates), end(cD@coordinates), ".", "*", ".", paste("ID=L", 1:nrow(cD), sep = "")), sep = "\t", append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE, file = file))
  }
