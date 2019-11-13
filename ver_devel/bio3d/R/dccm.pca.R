"dccm.pca" <-
  function(x, pc = NULL, method=c("pearson", "lmi"), ncore = NULL, ...) {
   if (missing(x) || !"pca" %in% class(x))
     stop("dccm.pca: must supply a 'pca' object, i.e. from 'pca.xyz'")

   method = match.arg(method)
    
   modes = pc

   ## Check for multiple cores
   ncore = setup.ncore(ncore)

   ## Set modes to be included
   if(is.null(modes))
      modes <- 1:length(x$L)
  
   ## If modes are negative, take modes complementary to them
   if( any(!is.numeric(modes)) || 
       any(!(abs(modes) %in% c(1:length(x$L)))) ||
       !(all(modes>0) || all(modes<0)) )
      stop("Incorrect mode index")
   if(all(modes < 0)) {
      modes <- setdiff(c(1:length(x$L)), abs(modes))
      if(length(modes) == 0)
         stop("No mode is selected")
   } 

   modes <- unique(modes)
   nmodes <- length(modes)
    
   ## Calc variance-covariance matrix over a subset of modes
   vcovmat <- function(r.inds, pca, vcov.mat = 0) {
     for ( i in seq_along(r.inds) ) {
       vcov.mat <- vcov.mat + (pca$U[, r.inds[i]] %o% pca$U[, r.inds[i]]) * pca$L[r.inds[i]]
       .update.pb(pb)
     }
     return(vcov.mat)
   }

   ## Calculate variance-covariance matrix first ## 
   ## If contain $z, straightforward
   if(!is.null(x$z)) {

      q = x$z[, modes] %*% t(x$U[, modes])
      vcov = cov(q)

   } else {

      ## Initialize progress bar
      pb <- .init.pb(ncore, min=0, max=nmodes)

      if(ncore > 1) {   # Parallel

         jobid <- rep(1:ncore, ceiling(nmodes/ncore))
         jobid <- jobid[1:nmodes]

         ltv <- mclapply(1:ncore, function(i) {
            j <- which(jobid %in% i)
            if(length(j) > 0) {
               m <- vcovmat(modes[j], x)
               m <- m[lower.tri(m, diag = TRUE)]
            } else {
               m = 0
            }
            return(m)
         } )
         
         ltv <- colSums(do.call(rbind, ltv))
         vcov <- matrix(0, nrow(x$U), nrow(x$U))
         vcov[lower.tri(vcov, diag = TRUE)] <- ltv
         vcov <- vcov + t(vcov)
         diag(vcov) <- diag(vcov) / 2

      } else {       # Serial

         vcov <- vcovmat(modes, x)

      }
      .close.pb(pb)
   }
   
   corr.mat <- .cov2dccm(vcov, method = method, ncore = ncore)
   return(corr.mat)
}
