"plot.enma" <-
  function(x, xlab="Residue Position", ylab="Fluctuations",
           col=NULL, cex=0.8, mar = c(6, 5, 2, 2),
           ...) {

    if(!inherits(x, "enma"))
      stop("provide a enma object as obtained from 'nma.pdbs'")
    
    opar <- par(no.readonly = TRUE)
    on.exit(par(opar))
    par(cex = cex, mar = mar)

    if(is.null(col))
      col <- seq(1, nrow(x$fluctuations))
    
    plot.bio3d(x$fluctuations[1,], xlab=xlab, ylab=ylab,
               ylim=c(0,max(x$fluctuations, na.rm=TRUE)),
               col=col[1], cex=cex, ... )
    
    for(i in 2:nrow(x$fluctuations)) {
      lines( x$fluctuations[i,], col=col[i] )
    }
    
  }
