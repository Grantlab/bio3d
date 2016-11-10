`get.seq` <-
function(ids, outfile="seqs.fasta", db="nr", verbose=FALSE) {
  ## Download FASTA format sequences from the NR or
  ## SWISSPROT/UNIPROT databases via their gi or
  ## SWISSPROT identifer number

  if( !(db %in% c("nr", "swissprot", "uniprot")) )
    stop("Option database should be one of nr, swissprot or uniprot")


  ids <- unique(ids)
  
    if(db=="nr") {
      get.files <- paste("https://www.ncbi.nlm.nih.gov/",
                         "sviewer/viewer.fcgi?db=protein&val=",
                         ids,"&report=fasta&retmode=text", sep="")

##    ## Old pre Oct-18th-2010 format URL
##    get.files <- paste("http://www.ncbi.nlm.nih.gov/entrez/",
##                       "viewer.fcgi?db=protein&val=",
##                       ids, "&dopt=fasta&sendto=t", sep="")

  } else {
#    if(any(nchar(ids) != 6)) {
#      warning("ids should be standard 6 character SWISSPROT/UNIPROT formart: trying first 6 char...")
#      ids <- substr(basename(ids),1,6)
#    }
    ids <- unique(ids)
    get.files <- file.path("http://www.uniprot.org/uniprot",
                           paste(ids, ".fasta", sep="") )
  }

  ## Remove existing file
  if(file.exists(outfile)) {
    warning(paste("Removing existing file:",outfile))
    unlink(outfile)
  }

  if(!verbose & length(get.files) > 1)
    pb <- txtProgressBar(min=0, max=length(get.files), style=3)
  
  retry <- 0
  k <- 1
  rtn <- rep(NA, length(ids))
  tmp.fasta <- tempfile()
  while(k <= length(ids)) {
    res <- tryCatch({
      download.file( get.files[k], tmp.fasta, mode="w", quiet=!verbose)
    }, error = function(e) {
      return(1)
    })
    
    if(res != 0 & retry < 10) {
      message(paste("\nDownload failed. Re-trying URL:\n ", get.files[k]))
      retry <- retry + 1
      rtn[k] <- res
      k <- k
      Sys.sleep(0.5)
    }
    else {
      retry <- 0
      rtn[k] <- res
      if(res == 0) {
        test <- try(aln <- read.fasta(tmp.fasta), silent=TRUE)
        if(inherits(test, "try-error")) {
           warning(paste("Failed to get sequence for", ids[k], "(check the ID)"))
           rtn[k] <- 1
        }
        else {
           write.fasta(aln, file=outfile, append=TRUE)
        } 
      }
      k <- k+1

      if(!verbose & length(get.files) > 1)
        setTxtProgressBar(pb, k)
    }
  }

  if(!verbose & length(get.files) > 1)
    cat("\n")
  
  rtn <- as.logical(rtn)
  names(rtn) <- ids
  if(all(!rtn)) {
    return(read.fasta(outfile))
  } else {
    warning("Not all downloads were sucesfull, see returned values")
    return(rtn)
  }
}


