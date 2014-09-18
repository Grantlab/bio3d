#  Compute the Cartesian coordinates of lattice vectors.

cell.coords <- function(...)
  UseMethod("cell.coords")

cell.coords.default <- function(x, digits = 3, ...)
{
  if(!is.numeric(x)) stop("'cell' must be numeric")
  if(length(x) != 6) stop("'cell' must be a vector of length 3")
  
  x[4:6] <- x[4:6]*pi/180
  
  M <- matrix(ncol=3,nrow=3)
  M[ ,1] <- c(x[1],0,0)
  M[ ,2] <- c(x[2]* cos(x[6]),x[2]*sin(x[6]),0)
  M[1,3] <-   x[3]* cos(x[5])
  M[2,3] <-   x[3]*(cos(x[4])-cos(x[5])*cos(x[6]))/sin(x[6])
  M[3,3] <-   x[3]*sqrt(1+2*cos(x[4])*cos(x[5])*cos(x[6])
                        -(cos(x[4]))^2-(cos(x[5]))^2-(cos(x[6]))^2)/sin(x[6])
  M <- round(M, digits)
  dimnames(M) <- list(c("x","y","z"), c("a","b","c"))
  
  return(M)
}
