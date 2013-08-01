"dccm.nma" <-
  function(nma, nmodes=NULL) {
    
    if (missing(nma))
      stop("dccm.nma: must supply a 'nma' object, i.e. from 'nma'")
    if(!"nma" %in% class(nma))
      stop("dccm.nma: must supply 'nma' object, i.e. from 'nma'")
    
    ## Inner product between all pairs of residues
    cross.inner.prod <- function(a, b) {
      mat <- apply(a, 1, "%*%", t(b))
      return(mat)
    }
    
    if(!is.null(nma$frequencies)) {
      freqs <- nma$frequencies
    }
    else {
      freqs <- nma$force.constants
    }
    
    if(is.null(nmodes))
      nmodes <- length(nma$L)
    else {
      nmodes <- nmodes + nma$triv.modes
      if(nmodes>length(nma$L)) {
        warning("'nmodes' larger than the number of modes")
        nmodes <- length(nma$L)
      }
    }
   
    ## Calculate residue-wise inner products
    ptm <- proc.time()
    corr.mat <- matrix(0, nma$natoms, nma$natoms)
       
    pb <- txtProgressBar(min=(nma$triv.modes+1), max=nmodes + nma$natoms, style=3)
    for ( i in (nma$triv.modes+1):nmodes )  {
      mode <- matrix(nma$U[,i], ncol=3, byrow=TRUE)
      corr.mat <- corr.mat + (cross.inner.prod(mode, mode) / (freqs[i]**2))
      setTxtProgressBar(pb, i)
    }
    corr.mat[lower.tri(corr.mat)] <- t(corr.mat)[lower.tri(corr.mat)]
 
    ## Basis for normalization
    a <- NULL
    k <- i ## for ProgressBar !
    for ( i in 1:nrow(corr.mat) )   {
      tmp <- 0
      for ( j in (nma$triv.modes+1):nmodes )   {
        mode <- matrix(nma$U[,j], ncol=3, byrow=TRUE)
        tmp <- tmp +
          as.numeric((mode[i,] %*% mode[i,]) / (freqs[j]**2))
      }
      a <- c(a, sqrt(tmp))
      k <- k+1
      setTxtProgressBar(pb, k)
    }
    close(pb)
    bn <- a%o%a
    
    ## Normalized correlation matrix
    corr.mat <- corr.mat / bn
    class(corr.mat) <- c("dccm", "matrix")
    
    t <- proc.time() - ptm
    cat(" Done in", t[[1]], "seconds.\n")
    
    return(corr.mat)
  }
