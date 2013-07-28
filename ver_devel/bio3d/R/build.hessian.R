"build.hessian" <-
  function(xyz, pfc.fun, normalize=TRUE)  {

    if (!is.function(pfc.fun))
      stop("build.hessian: 'pfc.fun' must be a function")
    
    ## Coordinates
    x <- xyz[ seq(1, length(xyz), by=3) ]
    y <- xyz[ seq(2, length(xyz), by=3) ]
    z <- xyz[ seq(3, length(xyz), by=3) ]
    natoms <- length(x)
    
    ## Distance matrix and pair force constants (ff dependent)
    dist.mat <- dm.xyz(xyz, mask.lower=F)
    pfc <- matrix(mapply(pfc.fun, dist.mat), ncol=natoms)

    ## Allocate a 3Nx3N matrix
    H <- matrix(0, nrow=3*natoms, ncol=3*natoms)

    ## Double for-loop - re-implement for effiency
    inds <- seq(1, natoms*3, by=3)
    for ( i in 1:natoms ) {
        m <- inds[i]

        for ( j in 1:natoms ) {
            n <- inds[j]
            k <- pfc[i,j]

            if (i!=j) {
                l <- -k
                a <- x[j]-x[i]
                b <- y[j]-y[i]
                c <- z[j]-z[i]

                ## normalize distance vector
                if (normalize)
                    v <- normalize.vector(c(a,b,c))
                else
                    v <- c(a,b,c)
            }
            else {
                l <- 0
                v <- c(0,0,0)
            }
            
            H[m  , n  ] <- v[1] * v[1] * l
            H[m  , n+1] <- v[1] * v[2] * l
            H[m  , n+2] <- v[1] * v[3] * l

            H[m+1, n  ] <- v[2] * v[1] * l
            H[m+1, n+1] <- v[2] * v[2] * l
            H[m+1, n+2] <- v[2] * v[3] * l

            H[m+2, n  ] <- v[3] * v[1] * l
            H[m+2, n+1] <- v[3] * v[2] * l
            H[m+2, n+2] <- v[3] * v[3] * l

            ## Same as:
            ##H[m:(m+2), n:(n+2)] <- (v%o%v) * l
        }
    }


    i <- 0; j <- 0;

    ## Calculate diagonal super-elements of H
    for ( i in 1:nrow(H) ) {
        m <- sum(H[i, inds  ]) # x
        n <- sum(H[i, inds+1]) # y
        o <- sum(H[i, inds+2]) # z

        j <- rep(inds, each=3)[i]
        H[i, j  ] <- m * (-1)
        H[i, j+1] <- n * (-1)
        H[i, j+2] <- o * (-1)
    }

    return(H)
}
