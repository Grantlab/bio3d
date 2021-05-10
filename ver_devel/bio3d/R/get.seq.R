`get.seq` <-
function(ids, outfile="seqs.fasta", db="nr", verbose=FALSE) {
  ## Download FASTA format sequences from the NCBI nr,
  ## SWISSPROT/UNIPROT, or RCSB PDB databases via their gi,
  ## SWISSPROT identifer number, or PDB ids.
  oops <- requireNamespace("httr", quietly = TRUE)
  if(!oops)
    stop("Please install the httr package from CRAN")

  db <- tolower(db)
  if( !(db %in% c("nr", "swissprot", "uniprot", "pdb")) )
    stop("Option database should be one of nr, swissprot/uniprot, or pdb")

  db <- switch(db, nr='nr', swissprot='uniprotkb', 
                   uniprot='uniprotkb', pdb='pdb')

  ids <- toupper(ids)
  ids <- unique(ids)

  if(db == "nr") {
    baseUrl <- 'https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?db=protein'
  } else {
    baseUrl <- 'https://www.ebi.ac.uk/Tools/dbfetch/dbfetch'
 
    # check if API works
    url <- paste(baseUrl, '/dbfetch.databases', sep='')
    resp <- httr::GET(url)
    if(httr::http_error(resp))
      stop('Access to EMBL-EBI server failed.')
  }
  
  # fetch sequences
  ## Remove existing file
  if(file.exists(outfile)) {
    warning(paste("Removing existing file:",outfile))
    unlink(outfile)
  }
  
  if(verbose) {
    cat("Fetching sequences from\n\t", baseUrl, "\n\nPlease wait", sep="")
  } else {
    cat("Fetching... Please wait")
  }

  ## do multiple requests if # of sequences > nmax
  if(grepl('ebi', baseUrl)) {
     nmax = 100
  } else {
     nmax = 500
  }
  errorCount=0
  checkInterval = 3
  checkInterval2 = 300 # wait longer if blocked by servers
  n <- floor( (length(ids)-1)/nmax) + 1
  for(i in 1:n) {
    if(i>1) Sys.sleep(10) # sleep 10s before sending another request
    i1 <- (i-1)*nmax+ 1
    i2 <- ifelse(i*nmax>length(ids), length(ids), i*nmax)
    ids1 <- ids[i1:i2]
    if(db=="nr") {
      url <- paste(baseUrl, "&val=", paste(ids1, collapse=','), 
                   "&report=fasta&retmode=text&page_size=",nmax, sep='')
    } else {
      url <- paste(baseUrl, db, paste(ids1, collapse=','), 'fasta?style=raw', sep='/')
    }
    resp <- try(httr::GET(url), silent=TRUE)
    if(inherits(resp, 'try-error')) {
      text <- 'LOST CONNECTION'
    } else {
      text <- httr::content(resp, 'text', encoding='utf-8')
    }
    retry=0
    while((grepl('^LOST CONNECTION', text) || httr::http_error(resp) || 
           grepl('^ERROR', text) ||
           grepl('Nothing has been found', text)) && retry<3) {
      retry = retry + 1
      if(grepl('^LOST CONNECTION', text)) {
        if(verbose) {
          cat('\nLost connection to the URL:\n\t', url, '\n\nRetry ', retry, 
            ' in ', as.integer(checkInterval2*retry/60), ' minutes...', sep='')
        } else {
          cat('\nLost connection. Retry ', retry, ' in ', 
            as.integer(checkInterval2*retry/60), 'min...', sep='') 
        }
        Sys.sleep(checkInterval2 * retry)
      } else { 
        if(verbose) {
          cat('\nFetching sequences failed from the URL:\n\t', url, '\n\nRetry ', 
            retry, '...', sep='')
        } else {
          cat('\nFailed. Retry ', retry, '...', sep='')
        }
        Sys.sleep(checkInterval)
      }
      resp <- try(httr::GET(url), silent=TRUE)
      if(inherits(resp, 'try-error')) {
         text <- 'LOST CONNECTION'
      } else {
         text <- httr::content(resp, 'text', encoding='utf-8')
      }
    }
    if(grepl('^LOST CONNECTION', text) || httr::http_error(resp) || 
       grepl('^ERROR', text) ||
       grepl('Nothing has been found', text)) {
      errorCount = errorCount + 1
      if(errorCount==n)
        stop('No sequence found. Check the ID(s).')
    } else {
       ## BUGFIX: remove the space between uniprot ID and acc. no.
       ## and so both will be read via read.fasta() for checking below.
       if(db == "uniprotkb") {
         text <- gsub("(>[^ ]*) ", "\\1|", text)
       }
       cat(text, file=outfile, append=TRUE)
    }
    cat('.')
  }
  cat(' Done.\n')
  
  # check if all sequences are downloaded successfully.
  seqs <- read.fasta(outfile)
  if(db=="nr") {
    myids <- strsplit(seqs$id, split='\\|')
    myids <- sapply(myids, function(x) {
      ii <- match('pdb', x)
      if(!is.na(ii)) {
        x[ii+1] <- paste(x[ii+1], x[ii+2], sep='_')
      }
      paste(x, collapse='|')
    })
    rtn <- sapply(ids, function(x) !any(grepl(x, myids)))
    
    if(length(seqs$id) == length(ids)) {
      if(any(rtn)) {
        warning("Returned sequence IDs may be different from query.")
      }
      rtn = FALSE
    }
  } 
  else {
    rtn <- sapply(ids, function(x) !any(grepl(x, seqs$id)))
  }
  
  if(all(!rtn)) {
    return(seqs)
  } else {
    warning("Not all downloads were successful. See returned values (TRUE=possibly failed).")
    return(rtn)
  }
}


