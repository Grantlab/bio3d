#' Calculate the distance between DCCMs and covariance matrices
#'
#' Compare two DCCMs/covariance matrices or, if provided with 
#' 3-D arrays, repeat comparison along the third dimension.
#'
#' @details
#' If \code{a} and \code{b} are matrices, return one numeric value.
#'
#' If \code{a} is an array and \code{b} is NULL, return a symmeric matrix 
#' containing distance between each other matrix of \code{a} along 
#' the third dimension of \code{a}.
#'
#' If \code{a} and \code{b} are arrays and \code{all.pairs=TRUE}, 
#' return an asymmeric matrix containing distance between each matrix 
#' of \code{a} and each matrix of \code{b}, along their corresponding
#' third dimension.
#' 
#' If \code{a} and \code{b} are arrays and \code{all.pairs=FALSE}, 
#' return a vector containing distance between each corresponding matrix 
#' of \code{a} and \code{b}. Dimensions of \code{a} and \code{b} must be
#' the same. 
#'
#' @param a a symmetric numeric matrix that is obtained from e.g. 
#'    \code{\link{dccm}} or \code{\link{cov}}, or a 3-D array containing
#'    a set of such matrices.
#' @param b the second matrix or array to compare with \code{a}. If NULL,
#'    comparison will be done within \code{a}. See details below.
#' @param all.pairs logical, if TRUE all pairwise distances between \code{a} 
#'        \code{b} along the third dimension are computed. If FALSE only the
#'         distances between corresponding matrices of \code{a} and \code{b}
#'         are computed. See details below.
#' @param ... additional arguments passed to \code{\link{dist}}.
#'
#' @return a numeric scalar, vector, or matrix depending on dimensions of 
#' input \code{a} and \code{b} and the value of \code{all.pairs}. 
#' See details above.
#'
#' @seealso \code{\link{dccm}}, \code{\link{cov}}, \code{\link{bhattacharyya}}
#'
#' @examples
#'
#' \dontrun{
#'
#' } 
dist.dccm <- function(a, b = NULL, all.pairs = TRUE, ...) {

   # check dimensions of inputs
   ndim.a = length(dim(a))
   ndim.b = length(dim(b))
   if(!is.null(b)) {
      if(ndim.a < 2 || ndim.a > 3)
         stop("'a' must be a matrix or 3-D array")
      if(ndim.b < 2 || ndim.b > 3)
         stop("'b' must be a matrix or 3-D array")
      if(all.pairs) {
         if(any(dim(a)[1L:2L] != dim(b)[1L:2L]))
            stop("'a' and 'b' must contain matrices having the same dimensions")
      }
      else {
         if(any(dim(a) != dim(b)))
            stop("'a' and 'b' must have the same dimensions when all.pairs=FALSE")
      }
   } 
   else {
      if(ndim.a != 3)
         stop("'a' must be a 3-D array when 'b' is absent")
   }
   if(ndim.a == 2) a <- array(a, dim=c(dim(a), 1))
   if(ndim.b == 2) b <- array(b, dim=c(dim(b), 1))
  
   names.a <- dimnames(a)[[3L]]
   names.b <- dimnames(b)[[3L]]

   # check symmetry and determine dccm or vcov, based on the first matrix of a
   a.chk = a[,,1]
   lower = a.chk[lower.tri(a.chk)]
   upper = t(a.chk)[lower.tri(a.chk)]
   if(!isTRUE(all.equal(lower, upper)))
      stop("Input(s) is not symmetric matrix")
   if(isTRUE(all.equal(diag(a.chk), rep(1, nrow(a.chk)))))
      bDCCM = TRUE
   else 
      bDCCM = FALSE

   # extract lower triangular parts
   lower = t( apply(a, 3, function(x) x[lower.tri(x, diag=!bDCCM)]) )
   if(!is.null(b)) {
      lower2 = t( apply(b, 3, function(x) x[lower.tri(x, diag=!bDCCM)]) )
      lower = rbind(lower, lower2)
   }

   # I guess it is still faster computing for all possible pairs 
   # using dist() than computing demanded pairs using inline codes
   full.mat = as.matrix(dist(lower, ...))
   
   if(nrow(full.mat) == 2) {
      out <- full.mat[1, 2]
   }
   else {
      if(is.null(b)) {
         out <- full.mat
         if(is.null(names.a)) 
            rownames(out) <- paste("a", 1:dim(a)[3L], sep="")
         else 
            rownames(out) <- names.a
         colnames(out) <- rownames(out)
      }
      else {
         if(all.pairs) {
            out <- full.mat[1:dim(a)[3L], (dim(a)[3L]+1):ncol(full.mat)]
            if(is.null(names.a))
               rownames(out) <- paste("a", 1:dim(a)[3L], sep="")
            else 
               rownames(out) <- names.a
            if(is.null(names.b))
               colnames(out) <- paste("b", 1:dim(b)[3L], sep="")
            else 
               colnames(out) <- names.b
         }
         else {
            inds <- cbind(1:dim(a)[3L], (dim(a)[3L]+1):ncol(full.mat))
            out <- full.mat[inds]
         }
      }
   }
   return(out) 
}
