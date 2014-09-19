filter.dccm <- function(x, xyz, fac = NULL, cutoff.cij = 0.4, 
          collapse = TRUE, extra.filter = NULL, ...) {

   # check cij format
   cij <- x
   if("all.dccm" %in% names(cij)) {
      cij <- cij$all.dccm
   } else if(is.list(cij)) {
      cij <- array(unlist(cij), dim = c(dim(cij[[1]]), length(cij)))
   } else if(is.matrix(cij)) {
      cij <- array(cij, dim = c(dim(cij), 1))
   } else if(!is.array(cij)) {
      stop("Input x should be an array/list containing correlation matrices")
   }
    
   # check factor vector for multiple networks construction
   if(!is.null(fac)) {
      if(!is.factor(fac)) fac = as.factor(fac)
   } else {
      fac = factor(rep("a", dim(cij)[3L])) 
   } 

   # check xyz for contact map calculation
   if(inherits(xyz, "pdbs")) {
      gaps.pos <- gap.inspect(xyz$xyz)
      xyz <- xyz$xyz[, gaps.pos$f.inds]
   }
   if(nrow(xyz) != dim(cij)[3L] && nlevels(fac) > 1)
      stop("xyz matrix doesn't match x. Set fac=NULL for single network construction")
   
   # convert cij to upper.tri matrix for internal use
   pcij <- apply(cij, 3, function(x) x[upper.tri(x)])
 
   ncij <- tapply(1:dim(cij)[3L], fac, function(i) {
      
      # contact map
      cm <- cmap(xyz[i, ], ...) 

      cij.min = apply(abs(pcij[, i]), 1, min)
      cij.max = apply(abs(pcij[, i]), 1, max)

      filter <- (cij.min >= cutoff.cij) | (cij.max >= cutoff.cij & cm[upper.tri(cm)]==1)
      
      if(!is.null(extra.filter))
         filter <- filter * extra.filter[upper.tri(extra.filter)] 
      
      ncij <- array(NA, dim=c(dim(cij[,,1]), length(i)))
      for(j in 1:dim(ncij)[3L]) {
         tcij <- cij[,,i[j]]
         tcij[upper.tri(tcij)] <- pcij[, i[j]] * filter
         tcij[lower.tri(tcij)] <- t(tcij)[lower.tri(tcij)]
         ncij[,,j] <- tcij
      }
      if(length(i) == 1) ncij <- ncij[,,1]
      return(ncij)
   } )
  
   if(collapse) ncij <- lapply(ncij, rowMeans, dims = 2) 
   if(nlevels(fac)==1) ncij <- ncij[[1]]

   return(ncij)
}
