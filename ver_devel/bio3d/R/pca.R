"pca" <- function(x, ...) {
  if(inherits(x, "matrix")) 
    class(x) <- c("matrix", "xyz")
  
  UseMethod("pca", x)
}
