findChunks <- function(alignments, gap, checkDuplication = TRUE)
{
  chunks <- Rle()
  maxChunk <- 0L

  for(chr in seqlevels(alignments)) {
    chraD <- which(seqnames(alignments) == chr)
    if(length(chraD) > 0)
      {
        chral <- ranges(alignments[chraD,])
        chunkNum <- c(0L, which(start(chral)[-1] - cummax(end(chral)[-length(chral)]) > gap))
        chunkID <- rep(1L:as.integer(length(chunkNum)) + maxChunk, diff(c(chunkNum, length(chral))))
        maxChunk <- as.integer(max(chunkID))
        chunks[chraD] <- chunkID
      }
  }
  
  values(alignments)$chunk <- chunks

  if(checkDuplication)
    {  
      chunkDup <- Rle()
      tags <- (as.character(values(alignments)$tag))
      rodal <- order(as.integer(chunks), tags)
      
      dups <- which(as.integer(chunks)[rodal] == c(as.integer(chunks)[rodal[-1]], Inf) & tags[rodal] == c(tags[rodal[-1]], Inf))
      
      chunkDups <- rep(FALSE, length(tags))
      chunkDups[c(dups, dups + 1L)] <- TRUE
      
      values(alignments)$chunkDup <- chunkDups[order(rodal)]
    }
      
  alignments
  }
