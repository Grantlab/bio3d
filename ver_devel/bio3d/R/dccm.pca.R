"dccm.pca" <-
  function(x, nmodes = NULL, ncore = NULL, ...) {
   if (missing(x) || !"pca" %in% class(x))
     stop("dccm.pca: must supply a 'pca' object, i.e. from 'pca.xyz'")

   ## Check for multiple cores
   ncore = setup.ncore(ncore)
   
   ## Set modes to be included
   if(is.null(nmodes)) {
     nmodes <- length(x$L)
   } else if(nmodes > length(x$L)) {
       warning("'nmodes' larger than the number of modes")
       nmodes <- length(x$L)
   }
   
   ## Calc variance-covariance matrix over a subset of modes
   vcovmat <- function(r.inds, pca, vcov.mat = 0) {
     for ( i in r.inds ) {
       vcov.mat <- vcov.mat + (pca$U[, i] %o% pca$U[, i]) * pca$L[i]
       if(ncore > 1) writeBin(1, fpb)
       else setTxtProgressBar(pb, i)
     }
     return(vcov.mat)
   }

   ## Calculate variance-covariance matrix first ## 
   ## If contain $z, straightforward
   if(!is.null(x$z)) {

      q = x$z[, 1:nmodes] %*% t(x$U[, 1:nmodes])
      vcov = cov(q)

   } else {

      ## Initialize progress bar
      pb <- txtProgressBar(min=1, max=nmodes, style=3)

      if(ncore > 1) {   # Parallel

         # For progress bar
         fpb <- fifo(tempfile(), open = "w+b", blocking = T)
         if(inherits(parallel:::mcfork(), "masterProcess")) {
            progress <- 0.0
            while(progress < nmodes && !isIncomplete(fpb)) {
               msg <- readBin(fpb, "double")
               progress <- progress + as.numeric(msg)
               setTxtProgressBar(pb, progress)
            }
            parallel:::mcexit()
         }
         ###################
 
         jobid <- rep(1:ncore, ceiling(nmodes/ncore))
         jobid <- jobid[1:nmodes]

         ltv <- mclapply(1:ncore, function(i) {
            j <- which(jobid %in% i)
            if(length(j) > 0) {
               m <- vcovmat(j, x)
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

         close(fpb)
         parallel:::readChildren()  # cleanup child process

      } else {       # Serial

         vcov <- vcovmat(1:nmodes, x)

      }
      close(pb)
   }
   
   corr.mat <- cov2dccm(vcov, ncore = ncore, ...)
   return(corr.mat)
}
