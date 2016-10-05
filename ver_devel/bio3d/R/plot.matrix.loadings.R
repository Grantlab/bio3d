#' Plot Residue-Residue Matrix Loadings 
#'
#' Plot residue-residue matrix loadings of a particular PC that is obtained from a 
#' principal component analysis (PCA) of cross-correlation or distance matrices.
#'
#' The function plots loadings (the eigenvectors) of PCA performed on a set of matrices 
#' such as distance matrices from an ensemble of crystallographic structures 
#' and residue-residue cross-correlations or covariance matrices derived from 
#' ensemble NMA or MD simulation replicates (See \code{\link{pca.array}} for detail). 
#' Loadings are displayed as a matrix with dimension the same as the input matrices 
#' of the PCA. Each element of loadings represents the proportion that the corresponding 
#' residue pair contributes to the variance in a particular PC. The plot can be used 
#' to identify key regions that best explain the variance of underlying matrices.
#'
#' @param x the results of PCA as obtained from \code{\link{pca.array}}. 
#' @param pc the principal component along which the loadings will be shown. 
#' @param resno numerical vector or \sQuote{pdb} object as obtained from \code{\link{read.pdb}}
#'    to show residue number on the x- and y-axis. 
#' @param sse a \sQuote{sse} object as obtained from \code{\link{dssp}} or \code{\link{stride}}, 
#'    or a \sQuote{pdb} object as obtained from \code{\link{read.pdb}} to show secondary
#'    structural elements along x- and y-axis.
#' @param mask.n the number of elements from the diagonal to be masked from output. 
#' @param ... additional arguments passed to \code{\link{plot.dccm}}.
#'
#' @return Plot and also returns a numeric matrix containing the loadings.
#'
#' @seealso 
#'      \code{\link{plot.dccm}}, \code{\link{pca.array}}
#'
#' @author Xin-Qiu Yao
#' 
#' @references 
#'    Skjaerven, L. et al. (2014) \emph{BMC Bioinformatics} \bold{15}, 399.
#'    Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
#' 
#' @examples
#' \dontrun{
#'    attach(transducin)
#'    gaps.res <- gap.inspect(pdbs$ali)
#'    sse <- bounds.sse(pdbs$sse[1, gaps.res$f.inds])
#'
#'    # calculate modes
#'    modes <- nma(pdbs, ncore=NULL)
#'
#'    # calculate cross-correlation matrices from the modes
#'    cijs <- dccm(modes, ncore=NULL)$all.dccm
#'
#'    # do PCA on cross-correlation matrices
#'    pc <- pca.array(cijs)
#'
#'    # plot loadings
#'    l <- plot.matrix.loadings(pc, sse=sse)
#'    l[1:10, 1:10]
#'
#'    # plot loadings with elements 10-residue separated from diagonal masked
#'    plot.matrix.loadings(pc, sse=sse, mask.n=10)
#'
#' }
plot.matrix.loadings <- function(x, pc=1, resno=NULL, sse=NULL, mask.n=0, ...) {
   
   if(!inherits(x, 'pca') && grepl('pca.array', x$call))
      stop('Input x must be a "pca" object obtained from "pca.array()".')

   args.plot.dccm <- formals(plot.dccm)
   dots <- list(...)
   args <- dots[names(dots) %in% names(args.plot.dccm)]
   if('segment.min' %in% names(dots))
      segment.min <- dots$segment.min
   else
      segment.min <- args.plot.dccm$segment.min
   if('show' %in% names(dots))
      show <- dots$show
   else
      show <- as.character(args.plot.dccm$show)[2]
   if(is.na(show)) show <- 'full'
   if(!'main' %in% names(args))
     args$main <- paste('Loadings of PC ', pc, 
       ' (', round(x$L[pc]/sum(x$L)*100, 1), '%)', sep='')

   ## compute the dim of matrix
   M <- nrow(x$U)
   N <- (1 + sqrt(1 + 8*M)) / 2
   is.wholenumber <-
      function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
   if(!is.wholenumber(N))
      stop('Wrong dimension of eigenvectors detected. Check the input x.')

   if(is.null(resno)) resno <- 1:N
   if(is.pdb(sse)) {
    sse <- unname(pdb2sse(sse))
   }
   else if(inherits(sse, 'sse')) {
    if(!is.null(sse$sse)) sse <- unname(sse$sse)
   }
   if(is.character(sse)) sse <- bounds.sse(sse)
   
   args$resno <- resno
   args$sse <- sse
   
   ## normalize
   x$U[, pc] <- x$U[, pc] / max(abs(x$U[, pc]))

   lmat <- matrix(0, N, N)
   lmat[upper.tri(lmat)] <- x$U[, pc]
   lmat[lower.tri(lmat)] <- t(lmat)[lower.tri(lmat)]
  
   ## mask diag 
   tmat <- matrix(1, N, N)
   tmat[diag.ind(tmat, n=mask.n)] <- 0
   tmat[lower.tri(tmat)] <- t(tmat)[lower.tri(tmat)]
   lmat <- lmat * as.vector(tmat)
  
   args$x <- lmat

   do.call(plot.dccm, args)

   ## add grids
   draw.sse.grid <- function(sse) {
      # vertical
      grid.segments( x0 = sse$start,
                     y0 = switch(show, full=1, upper=N, lower=1),
                     x1 = sse$start,
                     y1 = switch(show, full=N, upper=sse$start, lower=sse$start),
                     gp=gpar(col="gray80", lty=2, lwd=0.3), default.units = "native",
                     vp=vpPath("plot_01.toplevel.vp", "plot_01.panel.1.1.vp") )
      # horizental
      grid.segments( x0 = switch(show, full=1, upper=1, lower=N),
                     y0 = sse$start,
                     x1 = switch(show, full=N, upper=sse$start, lower=sse$start),
                     y1 = sse$start,
                     gp=gpar(col="gray80", lty=2, lwd=0.3), default.units = "native",
                     vp=vpPath("plot_01.toplevel.vp", "plot_01.panel.1.1.vp") )
   }

   if(!is.null(sse)) {
   	  if(length(sse$helix$start) > 0) {
        ## dont have a pdb$helix$length
        if( is.null(sse$helix$length) )
          sse$helix$length <- (sse$helix$end+1)-sse$helix$start

          inds <- which(sse$helix$length >= segment.min)
#          sse$helix$start <- match(sort(sse$helix$start[inds]), resno)
#          sse$helix$end <- match(sort(sse$helix$end[inds]), resno)
          draw.sse.grid(sse$helix)
      }
      if(length(sse$sheet$start) > 0) {
        ## dont have a pdb$sheet$length
        if( is.null(sse$sheet$length) )
          sse$sheet$length <- (sse$sheet$end+1)-sse$sheet$start

          inds <- which(sse$sheet$length >= segment.min)
#          sse$sheet$start <- match(sort(sse$sheet$start[inds]), resno)
#          sse$sheet$end <- match(sort(sse$sheet$end[inds]), resno)
          draw.sse.grid(sse$sheet)
      }
   }
   invisible( lmat )
}