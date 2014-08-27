"as.xyz" <- function(x) {
  if(is.vector(x))
    x = matrix(x, nrow=1)
  
  dims <- dim(x)
  if(!(dims[2L]%%3==0))
    stop("number of cartesian coordinates must be a multiple of 3")
  
  class(x) <- "xyz"
  return(x)
}

