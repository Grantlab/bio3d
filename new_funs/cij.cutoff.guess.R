# Estimate cij.cutoff as quantile Pr(cij<=cij.cutoff) = p
cij.cutoff.guess <- function(cij, p = NULL, ref.cutoff = NULL) {

   if(is.null(p))
      if(is.null(ref.cutoff)) p = seq(0.900, 0.995, 0.005)
      
   cijs <- cij
   if("all.dccm" %in% names(cijs)) cijs <- cijs$all.dccm
   if(is.array(cijs) && length(dim(cijs))==3)
      cijs <- do.call("c", apply(cijs, 3, list))
   if(!is.list(cijs))
      stop("cijs should be matrix, array(dim=3), or list")
   
   if(is.null(p)) {
      # return cumulative probability Pr(cij <= ref.cutoff)
      out <- sapply(cijs, function(x) 
          ecdf(abs(x[upper.tri(x)]))(ref.cutoff))
      out <- matrix(out, ncol=length(cijs))
      dimnames(out) <- list(c(paste("Pr(x<=", round(ref.cutoff, digits=2), ")", sep="")),
                            paste("Matrix", 1:length(cijs)))
   } else {
      # return quantile Pr(cij <= cij.cutoff) = p
      out <- sapply(cijs, function(x)
          quantile(abs(x[upper.tri(x)]), probs = p))
      out <- matrix(out, ncol=length(cijs))
      dimnames(out) <- list(paste("Cutoff for p=", round(p, digits=2), sep=""),
                            paste("Matrix", 1:length(cijs)))
   }
   return(out) 
}
