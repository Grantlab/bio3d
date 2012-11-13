
"overlap" <-
  function(pca, dv, num.modes=20) {
    
    if ( class(pca) == "pca" ) {
       ev <- pca$U
    } else {
       ev <- pca
    }

    if (nrow(ev)!=length(dv))
      stop("unequal vector lengths")

    if ( dim(ev)[2L] < num.modes ) {
       num.modes <- dim(ev)[2L]
       warning("num.modes larger than dimensions of pca")
    }
       
    normalize.vector <- function(v) {
      if(class(v)=='matrix')
	return( t(t(v)/sqrt( colSums(v**2)) ) )
      else	
	return( v/sqrt( sum(v**2) ) )
    }

    ev <- ev[,1:num.modes]
    overlap.values <- colSums(normalize.vector(ev) * normalize.vector(dv))**2
    ## t(normalize.vector(ev)) %*% normalize.vector(dv))**2 
        
    cum <- cumsum(overlap.values)
    out <- list(overlap=overlap.values, overlap.cum=cum) 
    
    return(out)
  }


