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
    baseUrl <- 'http://www.ebi.ac.uk/Tools/dbfetch/dbfetch'
 
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

  cat("Fetching... Please wait")
  
  ## do multiple requests if # of sequences > 200
  nmax = 200
  errorCount=0
  checkInterval = 3
  n <- floor( (length(ids)-1)/nmax) + 1
  for(i in 1:n) {
    i1 <- (i-1)*nmax+ 1
    i2 <- ifelse(i*nmax>length(ids), length(ids), i*nmax)
    ids1 <- ids[i1:i2]
    if(db=="nr") {
      url <- paste(baseUrl, "&val=", paste(ids1, collapse=','), 
                   "&report=fasta&retmode=text&page_size=",nmax, sep='')
    } else {
      url <- paste(baseUrl, db, paste(ids1, collapse=','), 'fasta', sep='/')
    }
    resp <- httr::GET(url)
    text <- httr::content(resp, 'text', encoding='utf-8')
    retry=0
    while((httr::http_error(resp) || grepl('^ERROR', text) ||
           grepl('Nothing has been found', text)) && retry<3) {
      retry = retry + 1
      if(verbose) cat('Fetching sequences failed. Retry ', retry, '...\n', sep='')
      Sys.sleep(checkInterval)
      resp <- httr::GET(url)
      text <- httr::content(resp, 'text', encoding='utf-8')
    }
    if(httr::http_error(resp) || grepl('^ERROR', text) ||
       grepl('Nothing has been found', text)) {
      errorCount = errorCount + 1
      if(errorCount==n)
        stop('No sequence found. Check the ID(s).')
    }
    cat(text, file=outfile, append=TRUE)
    cat('.')
  }
  cat(' Done.\n')
  
  # check if all sequences are downloaded successfully.
  seqs <- read.fasta(outfile)
  if(db=="nr") {
     if(length(seqs$id) == length(ids)) {
       rtn = FALSE
     } else {
       rtn <- sapply(ids, function(x) !any(grepl(x, seqs$id)))
     }
  } else {
    myids <- strsplit(seqs$id, split='\\|')
    myids <- sapply(myids, function(x) {
      ii <- match('pdb', x)
      if(!is.na(ii))
        x[ii+1] <- paste(x[ii+1], ifelse(is.na(x[ii+2]), 'A', x[ii+2]), sep='_')
      paste(x, collapse='|')
    })
    rtn <- sapply(ids, function(x) !any(grepl(x, myids)))
  }
  
  if(all(!rtn)) {
    return(seqs)
  } else {
    warning("Not all downloads were successful. See returned values.")
    return(rtn)
  }
}


