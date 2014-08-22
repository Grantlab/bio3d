"pca" <- function(x, ...) {
  if(inherits(x, "matrix")) 
    class(x) <- c("matrix", "xyz")
  
  UseMethod("pca", x)
}

"pca.pdbs" <- function(x, ...)
  UseMethod("pca", x)

"pca.3dalign" <- function(x, core.find=FALSE, ...) {
  ## Log the call
  cl <- match.call()
  
  if(core) {
    core <- core.find(x)
    x$xyz = pdbfit(x, core$c0.5A.xyz)
  }
  
  gaps.pos <- gap.inspect(x$xyz)
  pc <- pca.xyz(x$xyz[,gaps.pos$f.inds])

  pc$call=cl
  return(pc)
}
