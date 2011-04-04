findChunks <- function(aD, gap)
{
    fastUniques <- function(x)
      if(nrow(x) > 1) return(c(TRUE, rowSums(x[-1L,, drop = FALSE] == x[-nrow(x),,drop = FALSE]) != ncol(x))) else return(TRUE)
  
  aD@alignments$chunk <- NA
  maxChunk <- 0L
  for(chr in aD@chrs$chr) {
    chraD <- aD@alignments$chr == chr
    if(any(chraD))
      {
        chral <- subset(aD@alignments, subset = chraD, select = c("start", "end"))
        chunkNum <- c(0L, which(chral$start[-1] - cummax(chral$end[-nrow(chral)]) > gap))
        chunkID <- rep(1L:as.integer(length(chunkNum)) + maxChunk, diff(c(chunkNum, nrow(chral))))
        maxChunk <- as.integer(max(chunkID))
        aD@alignments$chunk[chraD] <- chunkID
      }
  }

  aD@alignments$chunkDup <- NA
  aligns <- subset(aD@alignments, select = c("chunk", "tag"))
  rodal <- order(aligns$chunk, as.factor(aligns$tag))
  aligns <- aligns[rodal,]
  chunkDups <- which(!fastUniques(aligns))
  chunkDups <- union(chunkDups, chunkDups - 1L)
  
  aligns <- data.frame(aligns, chunkDup = FALSE)
  aligns$chunkDup[chunkDups] <- TRUE
  aD@alignments$chunkDup[rodal] <- aligns$chunkDup
  
  fastUniques <- function(x)
    if(nrow(x) > 1) return(c(TRUE, rowSums(x[-1L,, drop = FALSE] == x[-nrow(x),,drop = FALSE]) != ncol(x))) else return(TRUE)
  
  aD
  }
