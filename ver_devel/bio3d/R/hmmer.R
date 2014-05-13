
"hmmer" <- function(seq, type='phmmer', db=NULL, verbose=TRUE, timeout=90) {
  oopsa <- require(XML)
  oopsb <- require(RCurl)
  if(!all(c(oopsa, oopsb)))
     stop("Please install the XML and RCurl package from CRAN")

 
  alnToSeq <- function(seq) {
    if(inherits(seq, "fasta"))
      seq <- seq$ali
    if(is.matrix(seq))
      seq <- as.vector(seq[1,!is.gap(seq[1,])])
    else
      seq <- as.vector(seq)
    return(paste(seq, collapse=""))
  }
  
  alnToStr <- function(seq) {
    if(!inherits(seq, "fasta"))
      stop("seq must be of type 'fasta'")
    tmpfile <- tempfile()
    write.fasta(seq, file=tmpfile)
    rawlines <- paste(readLines(tmpfile), collapse="\n")
    return(rawlines)
  }
  
  ##types.allowed <- c("phmmer", "hmmscan", "hmmsearch", "jackhmmer")
  types.allowed <- c("phmmer", "hmmscan")
  if(! type%in%types.allowed )
    stop(paste("Input type should be either of:", paste(types.allowed, collapse=", ")))
  
  if(type=="phmmer") {
    seq <- alnToSeq(seq)
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
  }
  if(type=="hmmscan") {
    seq <- alnToSeq(seq)
    if(is.null(db))
      db="pfam"
    db.allowed <- tolower(c("pfam", "gene3d", "superfamily", "tigrfam"))
    db <- tolower(db)
    if(!db%in%db.allowed)
      stop("db must be of type...")
    seqdb <- NULL
    hmmdb <- db
  }
  if(type=="hmmsearch") {
    ## requires stockholm input format
    seq <- alnToStr(seq)
    if(is.null(db))
      db="pdb"
    db.allowed <- tolower(c("pdb", "swissprot"))
    db <- tolower(db)
    if(!db%in%db.allowed)
      stop("db must be of type...")
    seqdb <- db
    hmmdb <- NULL
  }
  
  
  tmpurl <- paste('http://hmmer.janelia.org/search/', type, sep="")
  cat(tmpurl, "\n")
  curl.opts <- list(httpheader = "Expect:",
                    httpheader = "Accept:text/xml",
                    timeout = timeout, 
                    connecttimeout = timeout,
                    verbose = verbose,
                    followlocation = TRUE
                    )
  
  hmm <- postForm(tmpurl, hmmdb=hmmdb, seqdb=seqdb, seq=seq,
                  style = "POST",
                  .opts = curl.opts)
  
  xml <- xmlParse(hmm)
  data <- as.data.frame(t(xpathSApply(xml, '///hits', xpathSApply, '@*')),
                        stringsAsFactors=FALSE)
  return(data)
}
