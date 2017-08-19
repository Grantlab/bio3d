`plot.blast` <-
function(x, cutoff=NULL, cut.seed=NULL, cluster=TRUE, mar=c(2, 5, 1, 1), cex=1.5, ...) {

  ## b <- blast.pdb( pdbseq( read.pdb("4q21") ) )
  ## plot(b, 188)  ## cut.seed=110

    cl <- class(x)
    if("hit.tbl" %in% names(x))
        x <- x$hit.tbl


  panelplot <- function(z=x$mlog.evalue, ylab="-log(Evalue)", gp=gp, ...) {
    z=as.numeric(z)
    plot(z, xlab="", ylab=ylab, col=gps, ...)
    abline(v=gp, col="gray70", lty=3)

    pos=c(rep(4, length(gp))[-length(gp)],2)
    text(  gp, z[gp], 
          labels=paste0("Nhit=",gp ,", x=", round(z[gp])), 
          col="black", pos=pos, cex=cex, ...) ##"gray50"
  }


  ##- Setup plot arangment
  opar <- par(no.readonly = TRUE)
    on.exit(par(opar))

    if("blast" %in% cl)
        par(mfcol=c(4,1), mar=mar, cex.lab=cex)
    if("hmmer" %in% cl & "kg" %in% names(x))
        par(mfcol=c(4,1), mar=mar, cex.lab=cex)

  ##- Find the point pair with largest diff evalue
  dx <- abs(diff(x$mlog.evalue))
  dx.cut = which.max(dx)


  if(!is.null(cutoff)) {
    ##- Use suplied cutoff
    gps = rep(2, length(x$mlog.evalue))
    gps[ (x$mlog.evalue >= cutoff) ] = 1

  } else {

    if(cluster) {
      ## Ask USER whether to continue with clustering with many hits  
      nhit <- length(x$mlog.evalue)
      if(nhit > 1500) {
        cluster <- readline( paste0(" Note: ", nhit, 
          " hits, continue with TIME-CONSUMING clustering [y/n/q](n): ") )

        cluster <- switch(cluster, y=TRUE, yes=TRUE, q="QUIT", FALSE)
        if(cluster=="QUIT") { stop("user stop") }
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
      

  ##- Plot each alignment statistic with annotated grps
  panelplot(gp=gp.inds)
    panelplot(x$bitscore, ylab="Bitscore", gp=gp.inds)

    if("identity" %in% names(x))
        panelplot(x$identity, ylab="Identity", gp=gp.inds)
    if("alignmentlength" %in% names(x))
        panelplot(x$alignmentlength, ylab="Length", gp=gp.inds)
    if("kg" %in% names(x)) {
                
        tbl <- table(x$kg[inds], cut(x$score[inds], 20))
        tbl=tbl[, seq(ncol(tbl), 1), drop=FALSE]
        cols <- seq(1,nrow(tbl))
        barplot(tbl, col=cols, ylab="Frequency", border="grey50")
        box()
        
        legend("topleft", rownames(tbl), col=cols,  pch=15, ncol=3, 
               cex=cex*0.8, box.lwd = .5, box.lty=2, box.col = "grey50", bg = "white")

        
        tbl <- table(x$kg[!inds], cut(x$score[!inds], 20))
        tbl=tbl[, seq(ncol(tbl), 1), drop=FALSE]
        cols <- seq(1,nrow(tbl))
        barplot(tbl, col=cols, ylab="Frequency", border="grey50")
        
        legend("topleft", rownames(tbl), col=cols,  pch=15, ncol=3,
               cex=cex*0.8, box.lwd = .5, box.lty=2, box.col = "grey50", bg = "white")
        box()
    }


  ##- Return details of hits above cutoff
  out <- cbind("pdb.id"=x$pdb.id[inds], "acc"=x$acc[inds], "group"=gps[inds])
  rownames(out) <- which(inds)
  o <- list(hits=out, pdb.id=x$pdb.id[inds], acc=x$acc[inds], inds=inds)
  class(o) <- "blast" 
  return(invisible(o))
}
