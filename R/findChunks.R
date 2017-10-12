% modification on git from copied files
findChunks <- function(alignments, gap, checkDuplication = TRUE, justChunks = FALSE)
{
  chunks <- Rle(rep(NA, length(alignments)))
  chunkID <- maxChunk <- 0L

  for(chr in seqlevels(alignments)) {
    chraD <- which(seqnames(alignments) == chr)
    if(length(chraD) > 0) {
      if(length(chraD) == 1) {
        chunks[chraD] <- maxChunk + 1
      } else if(length(chraD) > 1)
        {
          chral <- ranges(alignments[chraD,])
          chunkNum <- c(0L, which(start(chral)[-1] - cummax(end(chral)[-length(chral)]) > gap))
          chunkID <- rep(1L:as.integer(length(chunkNum)) + maxChunk, diff(c(chunkNum, length(chral))))
          chunks[chraD] <- chunkID
        }
      maxChunk <- max(chunks, na.rm = TRUE)
    }
  }

  if(justChunks) return(chunks)
  
  values(alignments)$chunk <- chunks

  if(checkDuplication)
    {  
      chunkDup <- Rle()

      if("tag" %in% names(values(alignments))) {
        tags <- (as.character(values(alignments)$tag))
        rodal <- order(as.integer(chunks), tags)
        
        dups <- which(as.integer(chunks)[rodal] == c(as.integer(chunks)[rodal[-1]], Inf) & tags[rodal] == c(tags[rodal[-1]], Inf))      
        chunkDups <- rep(FALSE, length(tags))
        chunkDups[c(dups, dups + 1L)] <- TRUE      
        values(alignments)$chunkDup <- chunkDups[order(rodal)]
      } else values(alignments)$chunkDup <- FALSE
    }
      
  alignments
  }
