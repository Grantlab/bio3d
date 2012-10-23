

"rmsip" <-
  function(pca.a, pca.b, subset = 10, row.name=NULL, col.name=NULL) {
    
    if (nrow(pca.a$U)!=nrow(pca.b$U) )
      stop("unequal dimensions")
    
    if  (nrow(pca.a$U) < subset ) {
      subset <- nrow(pca.a$U)
      warning("subset larger than dimensions of pca")
    }
    
    normalize.vector <- function(v) {
      return( v/sqrt( colSums(v**2) ) )
    }
    
    o <- ( t(normalize.vector(pca.a$U[,1:subset]))
          %*% normalize.vector(pca.b$U[,1:subset]) )**2
    
    if (!is.null(row.name)) {
      rownames(o) <- paste(row.name, c(1:subset))
    }
    
    if (!is.null(col.name)) {
      colnames(o) <- paste(col.name, c(1:subset))
    }
    
    rmsip <- sqrt(sum(o)/subset)
    out <- list(overlap=round(o,3), rmsip=rmsip)
    
    return( out )
  }
