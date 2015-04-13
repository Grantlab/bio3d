phmmer_local <- function(seq) {
  exefile <- configuration$hmmer$exefile
  db <- configuration$hmmer$pdbseq
  
  status <- system(exefile, ignore.stderr = TRUE, ignore.stdout = TRUE)
  if(!(status %in% c(0,1))) {
    stop(paste("Launching external program failed\n",
               "  make sure '", exefile, "' is in your search path", sep=""))
  }
  
  if(!file.exists(db)) {
    stop("local database not found")
  }

  infile <- tempfile()
  outfile1 <- tempfile()
  outfile2 <- tempfile()  
  write.fasta(seq, file=infile)
  
  system(paste(exefile, "--noali --tblout", outfile2, "-o", outfile1, infile, db))
  
  return(.read_hmmerout(outfile2))
}


.read_hmmerout <- function(file) {

  trim <- function(s) {
    ##- Remove leading and trailing spaces from character strings
    s <- sub("^ +", "", s)
    s <- sub(" +$", "", s)
    s[(s=="")]<-""
    s
  }
  
  rawlines <- readLines(file)
  inds <- grep("#", rawlines)
  lines <- rawlines[-inds]
  
  acc <- trim(substr(lines, 1, 6))
  evalue <- as.numeric(trim(substr(lines, 65, 73)))
  score <- as.numeric(trim(substr(lines, 75, 80)))
  desc <- trim(substr(lines, 169, 220))
  return(data.frame(acc=acc, evalue=evalue, score=score, desc=desc,
                    stringsAsFactors=FALSE))
  
}

filter.hmmer <- function(x, cutoff=NULL, cut.seed=NULL, cluster=TRUE, value=FALSE) {
  
  if(is.null(x$evalue))
    stop("missing evalues")

  x$mlog.evalue=x$score

  ##- Find the point pair with largest diff evalue
  dx <- abs(diff(x$mlog.evalue))
  dx.cut = which.max(dx)
  
  
  if(!is.null(cutoff)) {
    ##- Use suplied cutoff
    gps = rep(2, length(x$mlog.evalue))
    gps[ (x$mlog.evalue >= cutoff) ] = 1

  } else {

    if(cluster) {
      ## avoid clustering with many hits  
      nhit <- length(x$mlog.evalue)
      if(nhit > 2000) {
        cluster <- FALSE
      }
    }

    if(is.null(cut.seed)) {
      ## Use mid-point of largest diff pair as seed for
      ##  cluster grps (typical PDB values are ~110)
      cut.seed = mean( x$mlog.evalue[dx.cut:(dx.cut+1)] )
    }

    if(cluster){
      ##- Partition into groups via clustering 
      ##  In future could use changepoint::cpt.var
      hc <- hclust( dist(x$mlog.evalue) )
      if(!is.null(cutoff)) { cut.seed=cutoff } 
      gps <- cutree(hc, h=cut.seed)
    } 

    if(!cluster || (length(unique(gps))==1)) {
      ##- Either we don't want to run hclust or hclust/cutree 
      ##   has returned only one grp so here we will divide   
      ##   into two grps at point of largest diff
      gps = rep(2, length(x$mlog.evalue))
      gps[1:dx.cut]=1
    }
  }

  gp.inds <- na.omit(rle2(gps)$inds)
  gp.nums <- x$mlog.evalue[gp.inds]


  
  cat("  * Possible cutoff values:   ", floor(gp.nums), "\n",
      "           Yielding Nhits:   ", gp.inds, "\n\n")

  if( is.null(cutoff) ) {
    ## Pick a cutoff close to cut.seed
    i <- which.min(abs(gp.nums - cut.seed))
    cutoff <- floor( gp.nums[ i ] )
  }

  inds <- x$mlog.evalue >= cutoff
  cat("  * Chosen cutoff value of:   ", cutoff, "\n",
      "           Yielding Nhits:   ", sum(inds), "\n")
  
  if(value) {
    return(x[inds, ])
  }
  else {
    out <- list(inds=inds, gp.inds=gp.inds, gps=gps, cutoff=cutoff)
    return(out)
  }
}

