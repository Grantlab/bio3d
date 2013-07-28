"normalize.vector" <-
  function(v, mass=NULL) {
    if(!is.null(mass)) {
      if (nrow(as.matrix(v)) != (length(mass)*3)) 
        stop("unequal vector lengths")
    }

    if(class(v)=='matrix')
      return(t( t(v) / sqrt(inner.prod(v,v,mass))))
    else
      return(v / sqrt(inner.prod(v,v,mass)))
  }
