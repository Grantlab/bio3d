seqbind <- function(..., blank = "-") {
  objs <- list(...)
  are.vec <- unlist(lapply(objs, is.vector))
  are.mat <- unlist(lapply(objs, is.matrix))
  if(!all(are.vec | are.mat))
    stop("'Can combine only vectors and/or matrices'")
  objs[are.vec] <- lapply(objs[are.vec], matrix, nrow = 1)
  max.col <- max(unlist(lapply(objs, ncol)))
  extend <- function(x, n, add)
    cbind(x, matrix(add, nrow=nrow(x), ncol=n-ncol(x)))  
  objs <- lapply(objs, extend, n = max.col, add = blank)
  return(do.call(rbind, objs))
}