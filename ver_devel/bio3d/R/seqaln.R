"seqaln" <-
function(aln, id=NULL, profile=NULL,
                   exefile = "muscle",
                   outfile = "aln.fa",
                   protein = TRUE,
                   seqgroup = FALSE,
                   refine = FALSE,
                   extra.args = "",
                   verbose = FALSE,
                   web.args = list()) {

  ## Log the call
  cl <- match.call()

  ## alignment to fasta object
  aln <- as.fasta(aln, id=id)

  ## nothing to align?
  if(!nrow(aln$ali) > 1 && is.null(profile)) {
    warning("nothing to align")
    aln$ali <- aln$ali[ , !is.gap(aln$ali), drop=FALSE]
    colnames(aln$ali) <- NULL
    return(aln)
  }
  
  if(!is.null(profile) & !inherits(profile, "fasta"))
    stop("profile must be of class 'fasta'")

  if(grepl("clustalo", tolower(exefile))) {
    prg <- "clustalo"
    ver <- "--version"
    
    if(!is.null(profile))
      args <- c("", "--profile1", "--in", "--out")
    else
      args <- c("--in", "--out")

    extra.args <- paste(extra.args,"--force")
    if(seqgroup)
      extra.args <- paste(extra.args, "--output-order=tree-order")
    else
      extra.args <- paste(extra.args, "--output-order=input-order")
    
    if(verbose)
      extra.args <- paste(extra.args,"--verbose")

    if(!is.null(profile) && length(grep("dealign", extra.args))==0)
      warning("profile alignment with clustalo: consider using extra.args='--dealign'")
      
    #if(protein)
    #  extra.args <- paste(extra.args,"--seqtype Protein")
    #else
    #  extra.args <- paste(extra.args,"--seqtype DNA")
  }
  else {
    prg <- "muscle"
    ver <- "-version"
    
    if(!is.null(profile))
      args <- c("-profile", "-in1", "-in2", "-out")
    else
      args <- c("-in", "-out")

    if(refine)
      extra.args <- paste(extra.args,"-refine")
    if(protein)
      extra.args <- paste(extra.args,"-seqtype protein")
    else
      extra.args <- paste(extra.args,"-seqtype dna")
  }
  
  ## Check if the program is executable
  os1 <- .Platform$OS.type
  status <- system(paste(exefile, ver),
                   ignore.stderr = TRUE, ignore.stdout = TRUE)

  if(!(status %in% c(0,1))) {
    if(length(web.args)==0) {
      stop(paste("You do not have ", prg, " installed/working locally on your machine.\n", 
        "  We can attempt to use the EBI webserver if you provide an email address (required by the EBI).\n",
        "  Please note that the EBI states (see their Terms of Use):\n",
        "     'Using fake e-mail address may result in your jobs being killed and your IP, Organisation or entire domain being black-listed.'\n", 
        sep=""))
    }
    cat('\n\nWill try to align sequences online...\n\n')
    
    default.web.args <- list(email='', title='', timeout=90)
    vnames <- names(web.args)
    vnames <- vnames[nzchar(vnames)]
    for(v in vnames) default.web.args[v] <- web.args[v]
    web.args <- default.web.args

    if(!is.null(profile))
      stop('Sequence-profile alignment is not supported by online web service.')
    if(!grepl('@', web.args$email))
      stop('A valid E-Mail address is required to use EMBL-EBI Web Service')
    
    naln <- .ebi_msa(aln, email=web.args$email, method=prg, protein=protein, 
                     title=web.args$title, timeout=web.args$timeout)
    if(!is.null(outfile)) write.fasta(naln, file=outfile)

  } else {

    ## Generate temporary files
    toaln <- tempfile()  
    write.fasta(aln, file=toaln)

    profilealn <- NULL
    if(!is.null(profile)) {
      profilealn <- tempfile()  
      write.fasta(profile, file=profilealn)
    }
    
    if(is.null(outfile))
      fa <- tempfile() 
    else
      fa <- outfile

    ## Build command to external program
    if(is.null(profile)) {
      cmd <- paste(exefile, args[1], toaln, args[2],
                   fa, extra.args, sep=" ")
    }
    else {
      cmd <- paste(exefile, args[1], args[2], profilealn, args[3], toaln, args[4],
                   fa, extra.args, sep=" ")
    }
    
    if(verbose)
      cat(paste("Running command:\n ", cmd , "\n"))
    
    ## Run command
    if (os1 == "windows")
      success <- shell(shQuote(cmd), ignore.stderr = !verbose, ignore.stdout = !verbose)
    else
      success <- system(cmd, ignore.stderr = !verbose, ignore.stdout = !verbose)
    
    if(success!=0)
      stop(paste("An error occurred while running command\n '",
                 exefile, "'", sep=""))
  
    naln <- read.fasta(fa, rm.dup=FALSE)

    ## Delete temporary files
    if(!is.null(profile))
      unlink(profilealn)
    unlink(toaln)
    if(is.null(outfile)) unlink(fa)
  }

  ## Re-group sequences to initial alignment order
  ## (muscle groups similar sequences by default)
  if(!seqgroup) {
    if(is.null(profile)) {
      ord <- match(aln$id, naln$id)
      naln$id <- naln$id[ord]
      naln$ali <- naln$ali[ord,]
    }
  }

  naln$call=cl
  return(naln)
}

## A client for EMBL-EBI multiple sequence alignment (MSA) Web Service.
## - 'x' must be a 'fasta' object.
## - 'email' must be a valid E-Mail address.
.ebi_msa <- function(x, email, method=c('muscle', 'clustalo'), protein=TRUE, title='', timeout=90) {
  method <- match.arg(method)

  oops <- requireNamespace("httr", quietly = TRUE)
  if(!oops)
    stop("Please install the httr package from CRAN")

#  oops <- requireNamespace("XML", quietly = TRUE)
#  if(!oops)
#    stop("Please install the XML package from CRAN")

  # Check number of sequences
  if(length(x$id)>500)
    stop("Number of sequences exceeds 500.") 

  baseUrl <- paste('http://www.ebi.ac.uk/Tools/services/rest/', method, sep='')

  # check if API works
  url <- paste(baseUrl, '/parameters', sep='')
  resp <- httr::GET(url)
  if(httr::http_error(resp))
    stop('EMBL-EBI msa-API request failed.')

  # submit a job
  url <- paste(baseUrl, '/run', sep='')
  seqs <- paste(paste('>', x$id, sep=''), apply(x$ali, 1, paste, collapse=''), 
                 sep='\n', collapse='\n')
  format <- switch(method, muscle='fasta', clustalo='fa')
  resp <- switch(method, 
    muscle = {
       httr::POST(url, body=list(email=email, title=title, format=format, 
             tree='none', order='aligned', sequence=seqs), encode='form')
       },
    clustalo = {
        httr::POST(url, body=list(email=email, title=title, 
          stype=ifelse(protein, 'protein', 'dna'),
          outfmt=format, dealign='true', sequence=seqs), encode='form')
    } )
  # check for errors
  if(httr::http_error(resp))
    stop('EMBL-EBI msa-API job submission failed.')

  # get job id 
  jobid <- httr::content(resp, 'text')
  cat('Job successfully submited (job ID: ', jobid, ')\n',
      'Waiting for job to finish...', sep='')
  
  # poll job status every 3s; stop if "error" obtained 3 times.
  # will also stopped if time out.
  url <- paste(baseUrl, '/status/', jobid, sep='')
  checkInterval = 3
  errorCount=0
  time = 0
  status <- 'PENDING'
  while((status %in% c('RUNNING', 'PENDING')) || 
        (status == 'ERROR' && errorCount < 2)) {
    status <- httr::content(httr::GET(url), 'text')
    if(status == 'ERROR') {
      errorCount <- errorCount + 1
    } else {
      if(errorCount > 0)
        errorCount <- errorCount - 1
    }
    if(status %in% c('RUNNING', 'PENDING', 'ERROR'))
      Sys.sleep(checkInterval)

    time <- time + checkInterval
    if(time >= timeout)
      stop(paste('\nConnection time out. Check your results from following URL:\n', 
        baseUrl, '/result/', jobid, '/aln-fasta', sep=''))
  }
  if(status != 'FINISHED') 
    stop('\nJob failed. Check your E-Mails for more information')

  cat('Done.\n')

  # check result type
  url <- paste(baseUrl, '/resulttypes/', jobid, sep='')
  resp <- httr::GET(url)
#  types <- XML::xmlToDataFrame(XML::xmlParse(httr::content(resp, 'text')), stringsAsFactors=FALSE)
  types <- httr::content(resp)
  types <- as.data.frame(t(sapply(types[[1]], function(x) x)), stringsAsFactors=FALSE)
  if(!any(grepl('fasta', types$identifier)))
    stop('Returned results do not contain FASTA format.')

  # get the results
  type <- types$identifier[grep('fasta', types$identifier)]
  url <- paste(baseUrl, '/result/', jobid, '/', type, sep='')
  resp <- httr::content(httr::GET(url), 'text')
  tfile <- tempfile()
  cat(resp, sep='', file=tfile)
  aln <- read.fasta(tfile)

  unlink(tfile)
  return(aln)
}
