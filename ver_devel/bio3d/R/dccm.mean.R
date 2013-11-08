dccm.mean <-
function(x, cutoff.sims = dim(x)[3], cutoff.cij = 0.4, ...) {

  ## Check input is a 3d array
  if (length(dim(x)) != 3) {
    stop("Input 'x' should be a NxNxS 3d array, where N is the number of atoms and S is the number of simulations")
  }
  ## Check input is built of simmetric matrices
  if (dim(x)[1] != dim(x)[2]) {
    stop("Input 'x' should be a NxNxS 3d array, where N is the number of atoms. The matrix uploaded is not NxNxS")
  }
  if (cutoff.sims > dim(x)[3]) {
    stop("The cutoff.sims number is greater than the number of simulations in the input matrix")
  }

  ## Filter by cutoff.cij and sum across simulations
  cut.cij.inds <- (abs(x) < cutoff.cij)
  count <- array(NA, dim = dim(x))
  count[!cut.cij.inds] = 1
  cij.sum <- apply(count, c(1:2), sum, na.rm = TRUE)

  ## Mask cij values below cutoff and average across simulations
  x[cut.cij.inds] = NA
  cij.ave <- apply(x, c(1:2), mean, na.rm = TRUE)

  ## Mask average values if below cutoff.sims
  cut.sims.inds <- (cij.sum < cutoff.sims)
  cij.ave[cut.sims.inds] = 0 ## Could use NA here
  
  class(cij.ave) = c("dccm", "matrix")
  return(cij.ave)
}
