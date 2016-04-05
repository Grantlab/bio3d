# Estimate cij.cutoff as quantile Pr(cij<=cij.cutoff) = p
cij.cutoff.guess <- function(cij, p = NULL, cmap = NULL, collapse = TRUE, collapse.method=c('max', 'median', 'mean')) {

   collapse.method <- match.arg(collapse.method)

   if(is.null(p)) p = seq(0.900, 0.995, 0.005)
      
   cijs <- cij
   if("all.dccm" %in% names(cijs)) cijs <- cijs$all.dccm
   if(is.array(cijs)) {
      if(length(dim(cijs))==3)
         cijs <- do.call("c", apply(cijs, 3, list))
      else
         cijs <- list(cijs)
   }
   if(!is.list(cijs))
      stop("cijs should be matrix, array(dim=3), or list")
   
   if(!is.null(cmap)) {
       cijs <- lapply(cijs, function(x) x*cmap)
   }

   # return quantile Pr(cij <= cij.cutoff) = p
   out <- sapply(cijs, function(x)
       quantile(abs(x[upper.tri(x)]), probs = p))
   out <- matrix(out, ncol=length(cijs))

   c0 <- seq(0, 1, 0.05)
   if(collapse) {
      out <- switch(collapse.method, 
        'mean' = 
#           sapply(rowMeans(out), function(x) c0[which.min(abs(x-c0))]) ,
           sapply(rowMeans(out), function(x) c0[sum(x>=c0)]),
        'median' = 
           sapply(apply(out, 1, median), function(x) c0[sum(x>=c0)]),
#           sapply(apply(out, 1, median), function(x) c0[which.min(abs(x-c0))]),
        'max' = 
           sapply(apply(out, 1, max), function(x) c0[sum(x>=c0)])
#           sapply(apply(out, 1, max), function(x) c0[which.min(abs(x-c0))])
      )
      names(out) <- paste("Cutoff (p=", round(p, digits=3), ")", sep="")
   } else {
      dimnames(out) <- list(paste("Cutoff (p=", round(p, digits=3), ")", sep=""),
                            paste("Matrix", 1:length(cijs)))
   }

   return(out) 
}
