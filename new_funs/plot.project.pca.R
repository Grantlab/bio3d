plot.project.pca <- function(pca, xyz, col = "gray", fit = FALSE, 
                           breaks = 60, ...) {

   require(grid)
   oops <- requireNamespace("lattice", quietly = TRUE)
   if(!oops)
      stop("Please install lattice package from CRAN")
   
   if(fit) {
      # Assume reference is the first structure in PCA
      q = as.vector(pca$U %*% pca$z[1, ]) + pca$mean
      xyz <- fit.xyz(q, xyz)
   }
   xyz <- project.pca(xyz, pca)[, 1:2]

   xr = range(c(pca$z[, 1], xyz[, 1]))
   yr = range(c(pca$z[, 2], xyz[, 2]))
 
   dots <- list(...)
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
      xlab = "PC1"
   }
   if("ylab" %in% names(dots)) {
      ylab = dots$ylab
      dots$ylab = NULL
   } else {
      ylab = "PC2"
   }
   
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
   data <- discretize2d(xyz[, 1], xyz[, 2], breaks, breaks)
   xx = seq(range(xyz[, 1])[1], range(xyz[, 1])[2], length.out = breaks + 1)
   yy = seq(range(xyz[, 2])[1], range(xyz[, 2])[2], length.out = breaks + 1)
   xx = xx[-length(xx)]
   yy = yy[-length(yy)]
   class(data) <- "matrix"

   # Plot trajectory sampling density
   p1 <- do.call(lattice::contourplot, c(list(x=data, region = TRUE, 
            labels = FALSE, contour = FALSE, xlab=xlab, ylab=ylab, 
            xlim = xlim, ylim = ylim, scales = list(alternating=1),
            row.values = xx, column.values = yy, aspect = 1,   
            panel = function(...) {
               lattice::panel.contourplot(...)
               lattice::panel.abline(h=0, v=0, lty=3, col="gray80")
            },
            at = seq(1, max(data), length.out = 255), colorkey = FALSE, 
            col.regions = colorRampPalette(blues9[-(1:3)])(255)), dots) )
   print(p1)

   # Add xray conformers 
   grid.points(pca$z[, 1], pca$z[, 2], pch=16, 
         gp = gpar(col = "black", cex = 0.6), 
         vp = vpPath("plot_01.toplevel.vp","plot_01.panel.1.1.vp") )
   grid.points(pca$z[, 1], pca$z[, 2], pch=16, 
         gp = gpar(col = col, cex=0.4), 
         vp = vpPath("plot_01.toplevel.vp","plot_01.panel.1.1.vp") )

   # Add start and end conformers of the trajectory
   grid.points(xyz[1, 1], xyz[1, 2], pch=16, 
         gp = gpar(col = "yellow", cex=0.4), 
         vp = vpPath("plot_01.toplevel.vp","plot_01.panel.1.1.vp") )
   grid.points(xyz[nrow(xyz), 1], xyz[nrow(xyz), 2], pch=16, 
         gp = gpar(col = "pink", cex=0.4), 
         vp = vpPath("plot_01.toplevel.vp","plot_01.panel.1.1.vp") )
  
}
