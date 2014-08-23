"pca.pdbs" <- function(x, core.find=FALSE, fit=FALSE, ...) {
  ## Log the call
  cl <- match.call()

  if(core.find & fit) {
    warning("incompatible arguments- neglecting 'fit=TRUE'")
    fit=FALSE
  }
  
  if(core.find) {
    core <- core.find(x)
    x$xyz = pdbfit(x, core$c0.5A.xyz)
  } else if(fit) {
     x$xyz = pdbfit(x)
  }
  
  gaps.pos <- gap.inspect(x$xyz)
  pc <- pca.xyz(x$xyz[,gaps.pos$f.inds], ...)

  pc$call=cl
  return(pc)
}
