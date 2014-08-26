"as.xyz" <- function(xyz) {
  if(is.vector(xyz))
    xyz = matrix(xyz, nrow=1)
  
  dims <- dim(xyz)
  if(!(dims[2L]%%3==0))
    stop("number of cartesian coordinates must be a multiple of 3")
  
  class(xyz) <- "xyz"
  return(xyz)
}

