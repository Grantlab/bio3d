"hmmer" <- function(seq, type='phmmer', db=NULL, verbose=TRUE, timeout=90) {
  cl <- match.call()
    
  oopsa <- require(XML)
  oopsb <- require(RCurl)
  if(!all(c(oopsa, oopsb)))
     stop("Please install the XML and RCurl package from CRAN")
 
  seqToStr <- function(seq) {
    if(inherits(seq, "fasta"))
      seq <- seq$ali
    if(is.matrix(seq)) {
      if(nrow(seq)>1)
        warning(paste("Alignment with multiple sequences detected. Using only the first sequence"))
      seq <- as.vector(seq[1,!is.gap(seq[1,])])
    }
    else
      seq <- as.vector(seq[!is.gap(seq)])
    return(paste(seq, collapse=""))
  }
  
  alnToStr <- function(seq) {
    if(!inherits(seq, "fasta"))
      stop("seq must be of type 'fasta'")
    tmpfile <- tempfile()
    write.fasta(seq, file=tmpfile)
    rawlines <- paste(readLines(tmpfile), collapse="\n")
    unlink(tmpfile)
    return(rawlines)
  }
  
  types.allowed <- c("phmmer", "hmmscan", "hmmsearch", "jackhmmer")
  if(! type%in%types.allowed )
    stop(paste("Input type should be either of:", paste(types.allowed, collapse=", ")))

  ## PHMMER (protein sequence vs protein sequence database)
  ## seq is a sequence
  if(type=="phmmer") {
    seq <- seqToStr(seq)
    if(is.null(db))
      db="pdb"
    db.allowed <- c("env_nr", "nr", "refseq", "pdb", "rp15", "rp35", "rp55",
                    "rp75", "swissprot", "unimes", "uniprotkb",
                    "uniprotrefprot", "pfamseq")
    db <- tolower(db)
    if(!db%in%db.allowed)
      stop("db must be of type...")
    
    seqdb <- db
    hmmdb <- NULL
    iter <- NULL
    rcurl <- TRUE
  }

  ## HMMSCAN (protein sequence vs profile-HMM database)
  ## seq is a sequence
  if(type=="hmmscan") {
    seq <- seqToStr(seq)
    if(is.null(db))
      db="pfam"
    db.allowed <- tolower(c("pfam", "gene3d", "superfamily", "tigrfam"))
    db <- tolower(db)
    if(!db%in%db.allowed)
      stop("db must be of type...")
    
    seqdb <- NULL
    hmmdb <- db
    iter <- NULL
    rcurl <- TRUE
  }

  ## HMMSEARCH (protein alignment/profile-HMM vs protein sequence database)
  ## seq is an alignment
  if(type=="hmmsearch") {
    if(!inherits(seq, "fasta"))
      stop("please provide 'seq' as a 'fasta' object")
    
    ##alnfile <- tempfile()
    ##seq <- write.fasta(seq, file=alnfile)
    seq <- alnToStr(seq)

    if(is.null(db))
      db="pdb"
    db.allowed <- tolower(c("pdb", "swissprot"))
    db <- tolower(db)
    if(!db%in%db.allowed)
      stop("db must be of type...")
    
    seqdb <- db
    hmmdb <- NULL
    iter <- NULL
    rcurl <- TRUE
  }

  ## JACKHMMER (iterative search vs protein sequence database)
  ## seq can be sequence, HMM, or multiple sequence alignment
  ## HMM not implemented here yet
  if(type=="jackhmmer") {
    if(!inherits(seq, "fasta"))
      stop("please provide 'seq' as a 'fasta' object")
    
    ##alnfile <- tempfile()
    ##seq <- write.fasta(seq, file=alnfile)
    seq <- alnToStr(seq)
    
    if(is.null(db))
      db="pdb"
    db.allowed <- tolower(c("pdb", "swissprot"))
    db.allowed <- c("env_nr", "nr", "refseq", "pdb", "rp15", "rp35", "rp55",
                    "rp75", "swissprot", "unimes", "uniprotkb",
                    "uniprotrefprot", "pfamseq")
    db <- tolower(db)
    if(!db%in%db.allowed)
      stop("db must be of type...")
    seqdb <- db
    hmmdb <- NULL
    iter <- NULL
    rcurl <- TRUE
  }
  
  
  ## Make the request to the HMMER website
  url <- paste('http://hmmer.janelia.org/search/', type, sep="")
  if(rcurl) {
    curl.opts <- list(httpheader = "Expect:",
                      httpheader = "Accept:text/xml",
                      ##timeout = timeout, 
                      ##connecttimeout = timeout,
                      verbose = verbose,
                      followlocation = TRUE
                      )
    
    hmm <- postForm(url, hmmdb=hmmdb, seqdb=seqdb, seq=seq, 
                    style = "POST",
                    .opts = curl.opts,
                    .contentEncodeFun=curlPercentEncode, .checkParams=TRUE )
  }
  else {
    ## temporary workaround for hmmsearch, and jackhmmer
    ## Check if the program is executable
    os1 <- .Platform$OS.type
    status <- system(paste("curl", "--version"),
                     ignore.stderr = TRUE, ignore.stdout = TRUE)
    
    if(!(status %in% c(0,1)))
      stop(paste("Launching external program failed\n",
                 "  make sure '", "curl", "' is in your search path", sep=""))

    outfile <- tempfile()
    cmd = paste("curl", " -L -H 'Expect:' -H 'Accept:text/xml'",
      " -F seqdb=", db, " -F seq='<", alnfile, "' ",  url,
      " > ", outfile, sep="")

    if(verbose)
      cat("Running command: ", cmd, "\n")

    if (os1 == "windows")
      success <- shell(shQuote(cmd), ignore.stderr = !verbose, ignore.stdout = !verbose)
    else
      success <- system(cmd, ignore.stderr = !verbose, ignore.stdout = !verbose)

    if(success!=0)
      stop(paste("An error occurred while running command\n '",
                 "curl", "'", sep=""))
    
    hmm <- readLines(outfile)
    unlink(outfile)
  }
  
  xml <- xmlParse(hmm)
  data <- as.data.frame(t(xpathSApply(xml, '///hits', xpathSApply, '@*')),
                        stringsAsFactors=FALSE)
  
  ## convert to numeric
  fieldsToNumeric <- c("evalue", "pvalue", "score", "archScore", "ndom", "nincluded",
                       "niseqs", "nregions", "nreported", "bias", "dcl", "hindex")
  inds <- which(names(data) %in% fieldsToNumeric)
  
  for(i in 1:length(inds)) {
    tryCatch({
      data[[inds[i]]] = as.numeric(data[[inds[i]]])
    },
    warning = function(w) {
      #print(w)
      return(data[[inds[i]]])
    },
    error = function(e) {
      #print(e)
      return(data[[inds[i]]])
    }
    )
  }
  
  class(data) <- c("hmmer", type, "data.frame")
  ##data$call <- cl
  return(data)
}


##".write.stockholm" <- function(aln, file=NULL)  {
##  if(!inherits(aln, "fasta"))
##    stop("aln must be of type 'fasta'")
##  
##  rawlines <- "# STOCKHOLM 1.0"
##  rawlines <- c(rawlines, paste("#=GF SQ", nrow(aln$ali)))
##
##  for( i in 1:nrow(aln$ali)) {
##    seq <- aln$ali[i,]
##    seq[is.gap(seq)]="-"
##    seq=paste(seq, collapse="")
##    tmpline <- paste(aln$id[i], seq, sep=" ")
##      
##    rawlines <- c(rawlines, tmpline)
##  }
##  rawlines <- c(rawlines, "//")
##  
##  if(is.null(file))
##    return(paste(rawlines, collapse="\n"))
##  else
##    write.table(rawlines, file=file, append=FALSE, quote=FALSE,
##                col.names=FALSE, row.names=FALSE)
##}
