`dccm` <- function(x, ...) {
  
  if(inherits(x, "matrix"))
    class(x) <- c("matrix", "xyz")
  else if(inherits(x, "array"))
    class(x) <- c("matrix", "mean")
  
  UseMethod("dccm")
}
