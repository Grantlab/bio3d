filter.dccm <- function(x, cutoff.cij = NULL, cmap = NULL, xyz = NULL, fac = NULL, 
                        cutoff.sims = NULL, collapse = TRUE, extra.filter = NULL, ...) {

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

   ## Check input is built of simmetric matrices
   if (dim(cij)[1] != dim(cij)[2]) {
     stop("Input 'x' should contain symmetric matrices.")
   }

   ## Check xyz and set cmap   
   if(is.null(cmap)) {
      if(is.null(xyz)) cmap = FALSE
      else cmap = TRUE
   }

   ## Check dimension if cmap is matrix
   if(is.matrix(cmap) && !all.equal(dim(cmap), dim(cij)[1L:2L]))
      stop("Input 'cmap' does not match x")

   ## Check cutoff.cij
   if(is.null(cutoff.cij)) {
      cutoff.cij <- .cij.cutoff.guess(cij, p = 0.95)
   }

   ## Inspect cij values with respect to cutoff.cij and contact map
   if(is.matrix(cmap) || isTRUE(cmap)) {

      # check factor vector for multiple networks construction
      if(!is.null(fac)) {
         if(!is.factor(fac)) fac = as.factor(fac)
      } else {
         fac = factor(rep("a", dim(cij)[3L])) 
      } 
  
      # check for calculating cmap 
      if(isTRUE(cmap)) {
         if(is.null(xyz))
            stop("xyz coordinates or a 'pdbs' object must be provided for contact map calculation")

         cmap.args <- list(...)

         if(inherits(xyz, "pdbs")) {
            gaps.pos <- gap.inspect(xyz$xyz)
            xyz <- xyz$xyz[, gaps.pos$f.inds]
            cmap.default <- list(dcut=10.0)
            cmap.args <- .arg.filter(cmap.default, cmap.xyz, ...)
         }
         if(nrow(xyz) != dim(cij)[3L] && nlevels(fac) > 1)
            stop(paste("Input 'xyz' doesn't match 'x'", 
              "\tSet fac=NULL for single network construction", sep='\n'))
      }

      # convert cij to upper.tri matrix for internal use
      pcij <- apply(cij, 3, function(x) x[upper.tri(x)])
 
      ncij <- tapply(1:dim(cij)[3L], fac, function(i) {
        
         # contact map
         if(isTRUE(cmap)) { 
            if(nlevels(fac) > 1)
                cm <- do.call("cmap.xyz", c(list(xyz=xyz[i, ]), cmap.args))
            else
                cm <- do.call("cmap.xyz", c(list(xyz=xyz), cmap.args)) 
         } else {
            cm <- cmap
         }
 
         cij.min = apply(abs(pcij[, i, drop=FALSE]), 1, min)
         cij.max = apply(abs(pcij[, i, drop=FALSE]), 1, max)
   
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
#         if(length(i) == 1) ncij <- ncij[,,1]
         return(ncij)
      } )
     
      if(collapse) ncij <- lapply(ncij, rowMeans, dims = 2) 
      if(nlevels(fac)==1) ncij <- ncij[[1]]
      if(is.matrix(ncij)) class(ncij) = c("dccm", "matrix")
      if(is.list(ncij)) {
         ncij <- lapply(ncij, function(x) {
            class(x) <- c("dccm", "matrix")
            x
         })
      }

      return(ncij)

   } else {

      # Filter cijs with cutoff.sims and return mean dccm (dccm.mean())

      if(is.null(cutoff.sims)) cutoff.sims = dim(cij)[3L]
   
      if (cutoff.sims > dim(cij)[3L] || cutoff.sims < 0) {
        stop("The cutoff.sims should be a number between 0 and N, where N is the the number of simulations in the input matrix")
      }

      ## Filter by cutoff.cij and sum across simulations
      cut.cij.inds <- (abs(cij) < cutoff.cij)
      count <- array(NA, dim = dim(cij))
      count[!cut.cij.inds] = 1
      cij.sum <- apply(count, c(1:2), sum, na.rm = TRUE)
    
      ## Mask cij values below cutoff and average across simulations
      cij[cut.cij.inds] = NA
      cij.ave <- apply(cij, c(1:2), mean, na.rm = TRUE)
    
      ## Mask average values if below cutoff.sims
      cut.sims.inds <- (cij.sum < cutoff.sims)
      cij.ave[cut.sims.inds] = 0 ## Could use NA here

      if(!is.null(extra.filter))
         cij.ave <- cij.ave * extra.filter
    
      class(cij.ave) = c("dccm", "matrix")
      return(cij.ave)
   }
}

# Estimate cij.cutoff as quantile Pr(cij<=cij.cutoff) = p
.cij.cutoff.guess <- function(cij, p = NULL, cmap = NULL, collapse = TRUE, collapse.method=c('max', 'median', 'mean')) {

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
