# MODELS
# A. simple
#    Remove edges with abs(cij) < cutoff.cij
# B. minimal
#    1. Remove edges with max(abs(cij)) < cutoff.cij
#    2. Remove edges with min(abs(cij)) < cutoff.cij && max(dij) > cutoff.dm
#    3. Remove edges with var(abs(cij)) >= cutoff.var
# C. full
#    1. Remove edges with max(abs(cij)) < cutoff.cij
#    2. Remove edges with min(abs(cij)) < cutoff.cij && max(dij) > cutoff.dm
# D. pca
#    1. Remove edges with max(abs(cij)) < cutoff.cij
#    2. Remove edges with min(abs(cij)) < cutoff.cij && max(dij) > cutoff.dm
#    3. Remove edges with var(abs(cij)) >= cutoff.var && intersect({pc-loading_n(ij) < cutoff.loading; n=1:nmodes})

cij.filter <- function(cij, inds = 1:dim(cij)[3L], xyz = NULL,
         model = c("simple", "minimal", "full", "pca", "dist"), nmodes = 3,
         cutoff.var = 0.0025, cutoff.cij = 0.4,
         cutoff.loading = 0.02, cutoff.dm = 15, cutoff.pcon = 1, extra.filter = NULL) {

   model <- match.arg(model)
   if(inherits(xyz, "3dalign")) {
      gaps.pos <- gap.inspect(xyz$xyz)
      xyz <- xyz$xyz[inds, gaps.pos$f.inds]
   }

   # check cij
   if("all.dccm" %in% names(cij)) cij <- cij$all.dccm
   if(is.list(cij)) cij <- array(unlist(cij), dim = c(dim(cij[[1]]), length(cij)))
   if(is.matrix(cij)) cij <- array(cij, dim = c(dim(cij), 1))
   if(any(is.na(inds))) inds <- 1:dim(cij)[3L]
  
   # convert cij to upper.tri matrix for internal use
   pcij <- apply(cij, 3, function(x) x[upper.tri(x)])
   pcij <- matrix(pcij[, inds], ncol=length(inds))
 
   filter <- switch(model,
     "simple" = abs(pcij) >= cutoff.cij,
     {
       if(is.null(xyz) || length(inds) == 1)
          stop("Non-simple model needs xyz and at least two cij matrices")

       # distance filter 
       dms <- apply(xyz, 1, function(x) {
          d <- dm.xyz(x)
          d[upper.tri(d)]} )
#       dm.min <- apply(dms, 1, min)
#       dm.max <- apply(dms, 1, max)
       pcon <- apply(dms, 1, function(x) sum(x<=cutoff.dm)/length(x))

       # cij variance
       cij.var <- apply(abs(pcij), 1, var)
       cij.min = apply(abs(pcij), 1, min)
       cij.max = apply(abs(pcij), 1, max) 

       fil <- switch(model, 
          "minimal" = {
            f <- cij.max < cutoff.cij
            f <- f | (cij.min < cutoff.cij & pcon < cutoff.pcon)
            f | cij.var >= cutoff.var
          },
          "full" = {
            f <- cij.max < cutoff.cij
            f | (cij.min < cutoff.cij & pcon < cutoff.pcon)
          },
          "pca" = {
            pca.cij <- pca.xyz(abs(t(pcij)), use.svd=TRUE)
            f <- cij.max < cutoff.cij
            f | (cij.min < cutoff.cij & pcon < cutoff.pcon)
            f2 <- TRUE
            for(i in 1:nmodes) f2 <- f2 & abs(pca.cij$U[, i]) < cutoff.loading
            f | (cij.var >= cutoff.var && f2)
         },
         "dist" = {
            pcon < cutoff.pcon
       } )
       matrix(rep(!fil, length(inds)), nrow=length(fil))
   } )
   
   if(!is.null(extra.filter))
      filter <- filter * extra.filter[upper.tri(extra.filter)]

   pcij <- pcij * filter
   new.cij <- array(NA, dim=c(dim(cij[,,1]), length(inds)))
   for(i in 1:length(inds)) {
      tcij <- cij[,,inds[i]]
      tcij[upper.tri(tcij)] <- pcij[, i]
      tcij[lower.tri(tcij)] <- t(tcij)[lower.tri(tcij)]
      new.cij[,,i] <- tcij
   }
   if(length(inds) == 1) new.cij <- new.cij[,,1]
   return(new.cij)
}
