`plot.pca` <-
function(x, pch=16, col=par("col"), cex=0.8, mar=c(4, 4, 1, 1), ...) {

  opar <- par(no.readonly=TRUE)
  on.exit(par(opar))

  par(mfrow=c(2, 2), cex=cex, mar=mar)
  par(pty="s")
  plot(x$z[,1],x$z[,2], type="p", pch=pch, xlab="PC1", ylab="PC2", col=col, ...)
  abline(h=0,col="gray",lty=2);abline(v=0,col="gray",lty=2)
  plot(x$z[,3],x$z[,2], type="p", pch=pch, xlab="PC3", ylab="PC2", col=col,...)
  abline(h=0,col="gray",lty=2);abline(v=0,col="gray",lty=2)
  plot(x$z[,1],x$z[,3], type="p", pch=pch, xlab="PC1", ylab="PC3", col=col,... )
  abline(h=0,col="gray",lty=2);abline(v=0,col="gray",lty=2)
  plot.pca.scree(x$L, ...)
}

