"build.hessian" <-
  function(xyz, pfc.fun, mass.weights=NULL, fc.weights=NULL, ncore=1 )  {

    if(missing(xyz))
      stop("build.hessian: 'xyz' coordinates must be provided")
    
    if (!is.function(pfc.fun))
      stop("build.hessian: 'pfc.fun' must be a function")

    if(!is.null(mass.weights)) {
      if((length(xyz)/3)!=length(mass.weights))
        stop("build.hessian: 'mass.weights' and number of atoms does not match")
    }
    
    ## Check for multicore package
    if(ncore>1) {
      oops <- require(multicore)
      if (!oops) {
        warning("multicore package missing")
        ncore <- 1
      }
    }

    build.submatrix <- function(xyz, rinds, inds.x, inds.y, inds.z,
                                mass.weights, natoms, fc.weights=NULL) {
      
      ## Sub-matrix of the full Hessian
      Hsm <- matrix(0, ncol=3*length(rinds), nrow=natoms*3)
      
      ## indices relating atoms and colums in the sub-hessian
      col.inds <- seq(1, ncol(Hsm), by=3)
      ## weight indices
      inds <- rep(1:natoms, each=3)
            
      for ( i in 1:length(rinds) ) {
        ## atom number
        n <- rinds[i]

        ## Calculate difference vectors and force constants
        diff.vect <- t(t(xyz) - xyz[n,])
        
        ##dists <- apply(diff.vect, 1, function(x) sqrt(sum(x**2)))
        dists <- sqrt(rowSums(diff.vect**2))  ## quicker !

        ## Previous version pfc.fun was not vectorized
        ##force.constants <- apply(diff.vect, 1, calc.fc, pfc.fun) * (-1)

        ## pfc.fun takes a vector of distances
        force.constants <- pfc.fun(dists)

        if(!is.null(fc.weights)) {
          force.constants <- force.constants * fc.weights[n,]
        }

        force.constants <- (-1) * force.constants / (dists**2)

        ## ensure that its not 'Inf'
        force.constants[n] <- 0

        ##diff.vect <- t(normalize.vector(t(diff.vect)))
        diff.vect[n,] <- 0 
        
        ## Hessian elements
        dxx <- diff.vect[,1] * diff.vect[,1] * force.constants
        dyy <- diff.vect[,2] * diff.vect[,2] * force.constants
        dzz <- diff.vect[,3] * diff.vect[,3] * force.constants
        
        dxy <- diff.vect[,1] * diff.vect[,2] * force.constants
        dxz <- diff.vect[,1] * diff.vect[,3] * force.constants
        dyz <- diff.vect[,2] * diff.vect[,3] * force.constants

        ## Place the elements
        m <- col.inds[i]
        
        ## Off-diagonals 
        Hsm[inds.x, m   ] <- dxx
        Hsm[inds.y, m+1 ] <- dyy
        Hsm[inds.z, m+2 ] <- dzz

        Hsm[inds.y, m   ] <- dxy
        Hsm[inds.z, m   ] <- dxz
        
        Hsm[inds.x, m+1 ] <- dxy
        Hsm[inds.z, m+1 ] <- dyz
        
        Hsm[inds.x, m+2 ] <- dxz
        Hsm[inds.y, m+2 ] <- dyz

        ## Diagonal super elements
        Hsm[inds.x[n], m] <- sum(Hsm[inds.x, m]) * (-1)
        Hsm[inds.y[n], m] <- sum(Hsm[inds.y, m]) * (-1)
        Hsm[inds.z[n], m] <- sum(Hsm[inds.z, m]) * (-1)

        Hsm[inds.x[n], m+1] <- sum(Hsm[inds.x, m+1]) * (-1)
        Hsm[inds.y[n], m+1] <- sum(Hsm[inds.y, m+1]) * (-1)
        Hsm[inds.z[n], m+1] <- sum(Hsm[inds.z, m+1]) * (-1)

        Hsm[inds.x[n], m+2] <- sum(Hsm[inds.x, m+2]) * (-1)
        Hsm[inds.y[n], m+2] <- sum(Hsm[inds.y, m+2]) * (-1)
        Hsm[inds.z[n], m+2] <- sum(Hsm[inds.z, m+2]) * (-1)

        ## Mass weight Hessian
        if(!is.null(mass.weights)) {
          m.tmp <- mass.weights[n] ## mass of atom n
          Hsm[,m:(m+2)] <- Hsm[,m:(m+2)] * (1/m.tmp)
          Hsm[,m:(m+2)] <- Hsm[,m:(m+2)] * (1/mass.weights[inds])
        }
      }
      return(Hsm)
    }

    ## Coordinates
    xyz <- matrix(xyz, ncol=3, byrow=TRUE)
    natoms <- nrow(xyz)

    if(!is.null(fc.weights)) {
      if(!is.matrix(fc.weights))
        stop("'fc.weights' must be a numeric matrix")
      
      if((nrow(fc.weights) != natoms) ||
         (ncol(fc.weights) != natoms) )
        stop("'fc.weights' must be numeric matrix with dimensions NxN")
    }
    
    ## Convenient indices for accessing the hessian
    inds.x <- seq(1, natoms*3, by=3)
    inds.y <- inds.x+1
    inds.z <- inds.x+2

    ## Divide the job to multiple cores
    ## Possibly mclapply instead of a for-loop is more efficient
    if(ncore>1) {
      thread.ids <- sort( rep(1:ncore, length.out=natoms) )
      ncore <- length(unique(thread.ids))
      jobs <- list()
    
      for ( i in 1:ncore ) {
        rinds <- which(thread.ids==i)
        jobs[[i]] <- multicore::parallel( build.submatrix(xyz, rinds,
                                     inds.x, inds.y, inds.z,
                                     mass.weights, natoms, fc.weights) )
      }
      
      res <- multicore::collect(jobs, wait=TRUE)
      
      H <- NULL
      for ( job in res ) {
        H <- cbind(H, job)
      }
    }
    else {
      rinds <- seq(1, natoms)
      H <- build.submatrix(xyz, rinds,
                           inds.x, inds.y, inds.z,
                           mass.weights, natoms, fc.weights)
    }

    return(H)
  }
