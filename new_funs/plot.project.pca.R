## col, color for crystal structures
plot.project.pca <- function(pca, xyz, pc.axes=c(1,2), col='gray', 
                           breaks=60, extra.xyz=NULL, 
                           col.extra.xyz='blue', ...) {

   require(grid)
   oops <- requireNamespace("lattice", quietly = TRUE)
   if(!oops)
      stop("Please install lattice package from CRAN")

   dots <- list(...)
   argnames.ppca <- names(formals(project.pca)) 
   args.ppca <- dots[names(dots) %in% argnames.ppca]
   dots <- dots[!names(dots) %in% argnames.ppca]

   xyz <- do.call(project.pca, c(list(data=xyz, pca=pca), args.ppca))

   ## check percentage of variance of MD captured by PC1 and PC2
   vc1 <- round(var(xyz[, 1])/sum(apply(xyz, 2, var))*100, 1)
   vc2 <- round(var(xyz[, 2])/sum(apply(xyz, 2, var))*100, 1)

   ## Variance of Xray captured by PCs
   vcx1 <- round(pca$L[pc.axes[1]]/sum(pca$L)*100, 1)
   vcx2 <- round(pca$L[pc.axes[2]]/sum(pca$L)*100, 1)

   xyz <- xyz[, pc.axes]
   z <- pca$z[, pc.axes]
   if(!is.null(extra.xyz)) {
     extra.xyz <- do.call(project.pca, c(list(data=extra.xyz, pca=pca), 
                        args.ppca))[, pc.axes]
     if(is.vector(extra.xyz)) extra.xyz <- t(extra.xyz)
     xr = range(c(z[, 1], xyz[, 1], extra.xyz[, 1]))
     yr = range(c(z[, 2], xyz[, 2], extra.xyz[, 2]))
   } else {
     xr = range(c(z[, 1], xyz[, 1]))
     yr = range(c(z[, 2], xyz[, 2]))
   }

   if("xlim" %in% names(dots)) {
      xlim = dots$xlim
      dots$xlim = NULL
   } else {
      xlim = c(xr[1] - diff(xr)*0.05, xr[2] + diff(xr)*0.05)
   }
   if("ylim" %in% names(dots)) {
      ylim = dots$ylim
      dots$ylim = NULL
   } else {
      ylim = c(yr[1] - diff(yr)*0.05, yr[2] + diff(yr)*0.05)
   }
   if("xlab" %in% names(dots)) {
      xlab = dots$xlab
      dots$xlab = NULL
   } else {
      xlab = paste('PC',pc.axes[1],' (Xray: ',vcx1,'%; MD: ',vc1, '%)', sep='')
   }
   if("ylab" %in% names(dots)) {
      ylab = dots$ylab
      dots$ylab = NULL
   } else {
      ylab = paste('PC',pc.axes[2],' (Xray: ',vcx2,'%; MD: ', vc2, '%)', sep='')
   }
   
#   if(fit) {
      # Assume reference is the 'mean' structure in PCA
#      q = as.vector(pca$U %*% pca$z[1, ]) + pca$mean
#      q = pca$mean
#      xyz <- fit.xyz(q, xyz)
#   }

   
   # From entropy package 
   discretize2d <- 
      function (x1, x2, numBins1, numBins2, r1 = range(x1), r2 = range(x2)) 
      {
          b1 = seq(from = r1[1], to = r1[2], length.out = numBins1 + 
              1)
          b2 = seq(from = r2[1], to = r2[2], length.out = numBins2 + 
              1)
          y2d = table(cut(x1, breaks = b1, include.lowest = TRUE), 
              cut(x2, breaks = b2, include.lowest = TRUE))
          return(y2d)
      }
   
   # Calculate density 
#   data <- discretize2d(xyz[, 1], xyz[, 2], breaks, breaks)
#   xx = seq(range(xyz[, 1])[1], range(xyz[, 1])[2], length.out = breaks + 1)
#   yy = seq(range(xyz[, 2])[1], range(xyz[, 2])[2], length.out = breaks + 1)
   data <- discretize2d(xyz[, 1], xyz[, 2], breaks, breaks, r1=xlim, r2=ylim)
   xx = seq(xlim[1], xlim[2], length.out = breaks + 1)
   yy = seq(ylim[1], ylim[2], length.out = breaks + 1)
   xx = xx[-length(xx)]
   yy = yy[-length(yy)]
   data <- data / sum(data) / ((diff(xlim)/breaks) * (diff(ylim)/breaks))
   class(data) <- "matrix"

   # Plot trajectory sampling density
   p1 <- do.call(lattice::contourplot, c(list(x=data, region = TRUE, 
            labels = FALSE, contour = FALSE, xlab=xlab, ylab=ylab, 
            xlim = xlim, ylim = ylim, scales = list(alternating=1, tck=c(1,0)), 
            row.values = xx, column.values = yy, aspect = 1,   
            panel = function(...) {
               lattice::panel.contourplot(...)
               lattice::panel.abline(h=0, v=0, lty=3, col="gray80")
            },
            at = seq(min(data[data>0]), max(data), length.out = 255), colorkey = FALSE, 
            col.regions = colorRampPalette(blues9[-c(1:2)])(255)), dots) )
   print(p1)

   # Add xray conformers 
   grid.points(z[, 1], z[, 2], pch=16, 
         gp = gpar(col = "black", cex = 0.6), 
         vp = vpPath("plot_01.toplevel.vp","plot_01.panel.1.1.vp") )
   grid.points(z[, 1], z[, 2], pch=16, 
         gp = gpar(col = col, cex=0.4), 
         vp = vpPath("plot_01.toplevel.vp","plot_01.panel.1.1.vp") )

   if(!is.null(extra.xyz)) {
     grid.points(extra.xyz[, 1], extra.xyz[, 2], pch=16, 
         gp = gpar(col = "black", cex = 0.6), 
         vp = vpPath("plot_01.toplevel.vp","plot_01.panel.1.1.vp") )
     grid.points(extra.xyz[, 1], extra.xyz[, 2], pch=16, 
         gp = gpar(col = col.extra.xyz, cex=0.4), 
         vp = vpPath("plot_01.toplevel.vp","plot_01.panel.1.1.vp") )
   }

   # Add start and end conformers of the trajectory
   grid.points(xyz[1, 1], xyz[1, 2], pch=17, 
         gp = gpar(col = "green", cex=0.4), 
         vp = vpPath("plot_01.toplevel.vp","plot_01.panel.1.1.vp") )
   grid.points(xyz[nrow(xyz), 1], xyz[nrow(xyz), 2], pch=17, 
         gp = gpar(col = "red", cex=0.4), 
         vp = vpPath("plot_01.toplevel.vp","plot_01.panel.1.1.vp") )
}
