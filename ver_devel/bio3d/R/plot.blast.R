`plot.blast` <-
function(x, cutoff=NULL, cut.seed=110, mar=c(4, 4, 1, 2), cex.lab=1.5, ...) {

  ## b <- blast.pdb( pdbseq( read.pdb("4q21") ) )
  ## plot(b, 188)
  
  nhit <- length(x$mlog.evalue)
  if(nhit > 2000) {
    continue <- readline(" More than 2000 hits, clustering may take some time, continue [y/n]:")
    continue <- ifelse(grepl("y",continue), TRUE, FALSE)
    if(!continue) { stop("user stop") }
  }

  ##- Setup Plot alignment stats overview
  opar <- par(no.readonly = TRUE)
  on.exit(par(opar))
  par(mfcol=c(4,1), mar=mar, cex.lab=cex.lab)

  
  
  ## Partition into groups
  hc <- hclust( dist(x$mlog.evalue) )
  if(!is.null(cutoff)) { cut.seed=cutoff } ## Need a better way to estimate grps
  gps <- cutree(hc, h=cut.seed)

  gp.mins <- rep(NA, max(gps))
  for(i in 1:max(gps)) {
    gp.mins[i] <- min(x$mlog.evalue[gps==i])
  }
  gp.end <- c(which(as.logical(diff(gps))), length(gps))
  gp.cut <- floor(gp.mins)
  gp.num <- cumsum( table(gps) )

  cat("  * Possible cutoff values include:\n\t\t", floor(gp.mins),
      "\n    Yielding Nhits:\n\t\t", gp.num, " \n\n")

  if( is.null(cutoff) ) {
    ## Pick a cutoff close to 110
    cut.ind <- which.min(abs(gp.cut - cut.seed))
    cutoff <- gp.cut[ cut.ind ]

    cat(" ** Chosen cutoff value of:\n\t\t", cutoff,
        "\n    Yielding Nhits:\n\t\t", gp.num[cut.ind], " \n")
      
  }
    
  plot(x$mlog.evalue, xlab="Hit No", ylab="-log(Evalue)", col=gps, ...)
  text(  gp.end, gp.mins, labels=gp.cut, col="gray50", pos=4, ...)
  ## abline(h=gp.cut, col="gray70", lty=3)
  abline(h=cutoff, col="red", lty=3)
  plot(x$bitscore, xlab="Hit No", ylab="Bitscore", col=gps, ...)
  plot(x$hit.tbl[,"identity"], xlab="Hit No", ylab="Identity", col=gps, ...)
  plot(x$hit.tbl[,"alignmentlength"], xlab="Hit No", ylab="Length", col=gps, ...)
  
  inds <- x$mlog.evalue >= cutoff
  out <- cbind("pdb.id"=x$pdb.id[inds], "gi.id"=x$gi.id[inds], "group"=gps[inds])
  rownames(out) <- which(inds)
  return(list(hits=out, pdb.id=x$pdb.id[inds], gi.id=x$gi.id[inds]))
  ##return(out)
}
