`blast.pdb` <-
function(seq, database="pdb") {
  if(inherits(seq, "fasta")) {
    if(is.matrix(seq$ali)) {
      if(nrow(seq$ali)>1)
        warning("Multiple sequences detected - using only the first sequence in object")
      seq <- as.vector(seq$ali[1,])
    }
    else {
      seq <- as.vector(seq$ali)
    }
  }
  
  ## Run NCBI blastp on a given 'seq' sequence against a given 'database'
  if(!is.vector(seq)) {
    stop("Input 'seq' should be a single sequence as a single or multi element character vector")
  }
  seq <- paste(seq, collapse="")
  
  if( !(database %in% c("pdb", "nr", "swissprot")) )
    stop("Option database should be one of pdb, nr or swissprot")

  ##- Submit
  urlput <- paste("http://www.ncbi.nlm.nih.gov/BLAST/Blast.cgi?CMD=Put&DATABASE=",
                  database,"&HITLIST_SIZE=20000&PROGRAM=blastp&CLIENT=web&QUERY=",
                  paste(seq,collapse=""),
                  sep="")


  txt <- scan(urlput, what="raw", sep="\n", quiet=TRUE)
  rid <- sub("^.*RID = " ,"",txt[ grep("RID =",txt) ])

  cat(paste(" Searching ... please wait (updates every 5 seconds) RID =",rid,"\n "))

  
  ##- Retrive results via RID code (note 'Sys.sleep()')
  urlget <- paste("http://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get",
                  "&FORMAT_OBJECT=Alignment",
                  "&ALIGNMENT_VIEW=Tabular",
                  "&RESULTS_FILE=on",
                  "&FORMAT_TYPE=CSV",
                  "&ALIGNMENTS=20000",
                  "&RID=",rid, sep="")


  raw  <- read.csv(urlget,
                   header = FALSE, sep = ",", quote="\"", dec=".",
                   fill = TRUE, comment.char="")

  
  ## Check for job completion (retrive html or cvs?)
  html <- 1
  while(length(html) == 1) {
    cat("."); Sys.sleep(5)
    
    raw  <- try(read.csv(urlget,
                     header = FALSE, sep = ",", quote="\"", dec=".",
                     fill = TRUE, comment.char=""), silent=TRUE)
    if(class(raw)=="try-error") { stop("No hits found: thus no output generated") }
    html <- grep("DOCTYPE", raw[1,])
  }

  
  colnames(raw) <- c("queryid", "subjectids", "identity", "positives",
                     "alignmentlength", "mismatches", "gapopens",
                     "q.start", "q.end", "s.start", "s.end",
                     "evalue", "bitscore")
  
  ## expand 'raw' for each hit in 'subjectids' (i.e. split on ";")
  rawm <- as.matrix(raw)
  
  eachsubject <- strsplit(rawm[,"subjectids"],";")
  subjectids  <- unlist(eachsubject)
  n.subjects  <- sapply(eachsubject, length)
  
  rawm <- apply(rawm, 2, rep, times=n.subjects)
  rawm[,"subjectids"] <- subjectids
  
  ## parse ids
  all.ids <- strsplit(subjectids, "\\|")
  gi.id  <- sapply(all.ids, '[', 2)
  pdb.id <- paste(sapply(all.ids, '[', 4),"_",sapply(all.ids, '[', 5),sep="")

  ## N.B. hack: zero evalues to arbirtrly high value!!
  mlog.evalue <- -log(as.numeric(rawm[,"evalue"]))
  mlog.evalue[is.infinite(mlog.evalue)] <- -log(1e-308)

  
  cat(paste("\n Reporting",length(pdb.id),"hits\n"))
  
  output <- list(bitscore=  as.numeric(rawm[,"bitscore"]),
                evalue =  as.numeric(rawm[,"evalue"]),
                mlog.evalue = mlog.evalue,
                gi.id = gi.id,
                pdb.id = pdb.id,
                hit.tbl = rawm,
                raw = raw)
  class(output) <- "blast"
  return(output)
}
