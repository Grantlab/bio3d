centres <- function(...)
  UseMethod("centres")

centres.default <- function(x, w = NULL, factor = NULL, unsplit = FALSE, ...){
  if(is.matrix(x))
    if(nrow(x)!=1)
      stop("'x' must be a single row 'xyz' matrix")
  x <- matrix(x, ncol=3, byrow=TRUE)
  if(is.null(w)) w <- rep(1, nrow(x))
  if(is.null(factor)) factor <- as.factor(rep(1, nrow(x)))
  if(length(w) != nrow(x))
    stop("length(w) must be equal to the number of atom inside 'x'")
  if(length(factor) != nrow(x))
    stop("length(factor) must be equal to the number of atom inside 'x'")
  
  w  <- split(w    , factor)
  x1 <- split(x[,1], factor)
  x2 <- split(x[,2], factor)
  x3 <- split(x[,3], factor)
  
  w.mean <- function(x, w)
    sum(x*w/sum(w, na.rm = TRUE))
  
  x1.mean <- mapply(w.mean, x1, w, SIMPLIFY = FALSE)
  x2.mean <- mapply(w.mean, x2, w, SIMPLIFY = FALSE)
  x3.mean <- mapply(w.mean, x3, w, SIMPLIFY = FALSE)
  
  if(unsplit){
    x1.mean <- unsplit(x1.mean, factor)
    x2.mean <- unsplit(x2.mean, factor)
    x3.mean <- unsplit(x3.mean, factor)
  }
  
  x1.mean <- unlist(x1.mean)
  x2.mean <- unlist(x2.mean)
  x3.mean <- unlist(x3.mean)

  x <- c(rbind(x1.mean, x2.mean, x3.mean))
  
  return(x)
}

centres.pdb <- function(x, w = atom2mass(x), factor = as.factor(x$atom[,"resno"]), unsplit = FALSE, ...){
  if(!is.pdb(x))
    stop("'x' must be an object of class 'pdb'")
  if(nrow(x$atom) != length(w))
    stop("length(w) must be equal to the number of atom inside 'x'")
  if(nrow(x$atom) != length(factor))
    stop("length(factor) must be equal to the number of atom inside 'x'")
  x <- centres.default(x$xyz, w, factor, unsplit)
  return(x)
}
