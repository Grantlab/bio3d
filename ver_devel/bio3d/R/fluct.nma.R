"fluct.nma" <-
  function(nma, mode.inds=NULL) {
    
    kb <- 0.00831447086363271
    pi <- 3.14159265359
    
    if(!"nma" %in% class(nma))
      stop("fluct.nma: must supply 'nma' object, i.e. from 'nma'")
    
    if("VibrationalModes" %in% class(nma))
      mass <- TRUE
    else
      mass <- FALSE
    
    if(is.null(mode.inds))
      mode.inds <- seq(nma$triv.modes+1, length(nma$L))
    
    f <- rep(0, nma$natoms)
    for ( i in mode.inds ) {
      mode <- matrix(nma$U[,i], ncol=3, byrow=TRUE)
      l <- apply(mode, 1, function(x) x%*%x)
      
      if(mass)
        l <- l/(nma$frequencies[i]**2)
      else
        l <- l/(nma$force.constants[i])
      f <- f+l
    }
    
    if(mass) {
      f <- f / nma$mass**2
      s <- 1/(2*pi)**2
      if(!is.null(nma$temp))
        s <- s*kb*nma$temp
      f <- f*s
    }
    else {
      if(!is.null(nma$temp))
        f <- f*kb*nma$temp
    }
    
    return(f)
  }
