seqbind <- function(..., blank = "-") {
  objs <- list(...)
  are.null <- unlist(lapply(objs, is.null))
  objs <- objs[!are.null]

  if(length(objs)==0)
    stop("Incompatible input")
  
  if(length(objs)==1)
    return(unlist(objs))
  
  are.vec <- unlist(lapply(objs, is.vector))
  are.mat <- unlist(lapply(objs, is.matrix))
  
  if(!all(are.vec | are.mat ))
    stop("'Can combine only vectors and/or matrices'")
  
  objs[are.vec] <- lapply(objs[are.vec], matrix, nrow = 1)
  max.col <- max(unlist(lapply(objs, ncol)))
  extend <- function(x, n, add)
    cbind(x, matrix(add, nrow=nrow(x), ncol=n-ncol(x)))  
  objs <- lapply(objs, extend, n = max.col, add = blank)
  return(do.call(rbind, objs))
}
