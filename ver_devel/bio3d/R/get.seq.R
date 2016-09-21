`get.seq` <-
function(ids, outfile="seqs.fasta", db="refseq", verbose=FALSE) {
  ## Download FASTA format sequences from the NCBI RefSeq,
  ## SWISSPROT/UNIPROT, or RCSB PDB databases via their gi,
  ## SWISSPROT identifer number, or PDB identifiers.
  oops <- requireNamespace("httr", quietly = TRUE)
  if(!oops)
    stop("Please install the httr package from CRAN")

  db <- tolower(db)
  if( !(db %in% c("refseq", "swissprot", "uniprot", "pdb")) )
    stop("Option database should be one of refseq, swissprot/uniprot, or pdb")
  db <- switch(db, refseq='refseqp', swissprot='uniprotkb', 
                   uniprot='uniprotkb', pdb='pdb')

  ids <- toupper(ids)
  ids <- unique(ids)

  baseUrl <- 'http://www.ebi.ac.uk/Tools/dbfetch/dbfetch'

  # check if API works
  url <- paste(baseUrl, '/dbfetch.databases', sep='')
  resp <- httr::GET(url)
  if(httr::http_error(resp))
    stop('EMBL-EBI server request failed.')

  # fetch sequences
  ## Remove existing file
  if(file.exists(outfile)) {
    warning(paste("Removing existing file:",outfile))
    unlink(outfile)
  }

  ## do multiple requests if # of sequences > 500
  errorCount=0
  checkInterval = 3
  n <- floor( (length(ids)-1)/500 ) + 1
  for(i in 1:n) {
    i1 <- (i-1)*500 + 1
    i2 <- ifelse(i*500>length(ids), length(ids), i*500)
    ids1 <- ids[i1:i2]
    url <- paste(baseUrl, db, paste(ids1, collapse=','), 'fasta', sep='/')
    resp <- httr::GET(url)
    text <- httr::content(resp, 'text')
    retry=0
    while((httr::http_error(resp) || grepl('^ERROR', text)) 
       && retry<3) {
      retry = retry + 1
      if(verbose) cat('Fetching sequences failed. Retry ', retry, '...\n', sep='')
      Sys.sleep(checkInterval)
      resp <- httr::GET(url)
      text <- httr::content(resp, 'text')
    }
    if(httr::http_error(resp) || grepl('^ERROR', text)) {
      errorCount = errorCount + 1
      if(errorCount==n)
        stop('No sequence found. Check the ID(s).')
    }
    cat(text, file=outfile, append=TRUE)
  }

  # check if all sequences are downloaded successfully.
  seqs <- read.fasta(outfile)
  myids <- strsplit(seqs$id, split='\\|')
  myids <- sapply(myids, function(x) {
    ii <- match('pdb', x)
    if(!is.na(ii))
      x[ii+1] <- paste(x[ii+1], ifelse(is.na(x[ii+2]), 'A', x[ii+2]), sep='_')
    paste(x, collapse='|')
    })
  rtn <- sapply(ids, function(x) !any(grepl(x, myids)))

  if(all(!rtn)) {
    return(seqs)
  } else {
    warning("Not all downloads were successful. See returned values.")
    return(rtn)
  }
}


