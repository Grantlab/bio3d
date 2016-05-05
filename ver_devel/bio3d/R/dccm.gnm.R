#' Dynamic Cross-Correlation from Gaussian Network Model
#'
#' Calculate the cross-correlation matrix from Gaussian network model normal 
#' modes analysis.
#'
#' This function calculates the cross-correlation matrix from Gaussian network
#' model (GNM) normal modes analysis (NMA) obtained from \code{gnm}. It returns 
#' a matrix of residue-wise cross-correlations whose elements, Cij, may be 
#' displayed in a graphical representation frequently termed a dynamical 
#' cross-correlation map, or DCCM. (See more details in \code{help(dccm.nma)}).
#'
#' @param x an object of class \sQuote{gnm} or \sQuote{egnm} as obtained from 
#'   \code{\link{gnm}}.
#' @param ... additional arguments (currently ignored).
#'
#' @return Returns a cross-correlation matrix.
#' 
#' @seealso \code{\link{gnm}}, \code{\link{dccm.nma}}, \code{\link{dccm.enma}},
#'   \code{\link{plot.dccm}}.
#'
#' @author Xin-Qiu Yao & Lars Skjaerven
#' 
#' @references 
#'    Bahar, I. et al. (1997) \emph{Folding Des.} \bold{2}, 173. 
#' 
#' @examples
#' ## Fetch stucture
#' pdb <- read.pdb( system.file("examples/1hel.pdb", package="bio3d") )
#'    
#' ## Calculate normal modes
#' modes <- gnm(pdb)
#'    
#' ## Calculate correlation matrix
#' cm <- dccm(modes)
#'
#' ## Plot correlation map
#' plot(cm, sse = pdb, contour = FALSE, col.regions = bwr.colors(20),
#'      at = seq(-1, 1, 0.1))
#'
#' @keywords analysis
"dccm.gnm" <- function(x, ...) {
  if (missing(x) || !inherits(x, 'gnm'))
    stop("dccm.gnm: must supply a 'gnm' object, i.e. from 'gnm.pdb()'")

  # variance-covariance matrix
  vcov <- cov.nma(x)
  
  # DCCM
  corr.mat <- vcov * 1/( sqrt(diag(vcov)) %*% t(sqrt(diag(vcov))) )
  class(corr.mat) <- c("matrix", "dccm")

  return(corr.mat)
}

#' @rdname dccm.gnm
"dccm.egnm" <- function(x, ...) {

  if (missing(x) || !inherits(x, 'egnm'))
  stop("dccm.egnm: must supply a 'egnm' object, i.e. from 'gnm.pdbs()'")
 
  # variance-covariance matrix
  vcovs <- cov.enma(x)

  # DCCM
  all.dccm <- apply(vcovs, 3, function(vcov) 
     vcov * 1/( sqrt(diag(vcov)) %*% t(sqrt(diag(vcov))) ) )
  all.dccm <- array(all.dccm, dim=dim(vcovs))

  avg.dccm <- rowMeans(all.dccm, dims=2)
  class(avg.dccm) <- c("matrix", "dccm")

  out <- list(all.dccm = all.dccm, avg.dccm = avg.dccm)
  return( out )

}
