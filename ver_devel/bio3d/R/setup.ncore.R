setup.ncore <- function(ncore) {
  if(is.null(ncore) || ncore > 1) {
    oops <- require(parallel)
    if(!oops) {
      if(is.null(ncore))
        ncore <- 1
      else
        stop("Please install the parallel package from CRAN\n\tor\n\tset ncore=1 for serial computation")
    }
    if(is.null(ncore))
      ncore = parallel:::detectCores()
    options(mc.cores = ncore)
  }
  return(ncore)
}
