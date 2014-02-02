setup.ncore <- function(ncore, bigmem = FALSE) {
  if(is.null(ncore) || ncore > 1) {
    oops <- require(parallel)
    if(!oops) {
      if(is.null(ncore)) {
        ncore <- 1
        bigmem <- FALSE
      } else {
        stop("Please install the parallel package from CRAN\n\tor\n\tset ncore=1 for serial computation")
      }
    }
    if(bigmem) {
      oops <- require(bigmemory)
      if(!oops) {
        if(is.null(ncore))
          ncore <- 1
        else
          stop("Please install the bigmemory package from CRAN for running with multicore")
      }
    }
    if(is.null(ncore))
      ncore = parallel:::detectCores()
    options(mc.cores = ncore)
  }
  return(ncore)
}
