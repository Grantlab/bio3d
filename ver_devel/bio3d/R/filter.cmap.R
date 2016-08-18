filter.cmap <- function(cm, cutoff.sims = NULL) {

  ## Convert list to array  
  if(is.list(cm)) {
      dims <- dim(cm[[1]])
      if(length(dims) != 2)
          stop("Input 'cm' should be a NxNxS 3d array or a list of NxN matrices,
      where N is the number of atoms and S is the number of simulations")
      
      if(length(unique(c(sapply(cm, dim)))) != 1)
          stop("Matrices in provided list have unequal dimensions")
      
      cm <- array(unlist(cm), dim = c(nrow(cm[[1]]), ncol(cm[[1]]), length(cm)))
  }
    
  ## Check input
  dims <- dim(cm)
  if(length(dims) != 3) {
      stop("Input 'cm' should be a NxNxS 3d array or a list of NxN matrices,
      where N is the number of atoms and S is the number of simulations")
  }

  if(dims[3L] < 2)
      stop("2 or more contact maps required")
  
  if(is.null(cutoff.sims))
      cutoff.sims <- dims[3L]

    if(!is.numeric(cutoff.sims) | cutoff.sims > dims[3L] | cutoff.sims < 1) {
        stop(paste("Input 'cutoff.sims' should be a number between 1 and", dims[3L], "\n  ",
                   "(i.e. the number of simulations upon which filtering is based)"))
    }

  ## Sum across simulations and filter by cutoff.sims de
  cm.sum <- apply(cm, c(1:2), sum, na.rm = TRUE)

  return(array(as.numeric(cm.sum >= cutoff.sims), dim = dim(cm.sum)))
}
