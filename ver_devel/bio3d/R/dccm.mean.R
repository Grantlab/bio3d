dccm.mean <-
function(cij, cutoff.sims = dim(cij)[3], cutoff.cij = 0.4) {

  ## Check if bio3d package is uploaded
  oops <- require(bio3d)
  if (!oops) {
    warning("bio3d package missing: Please install, see: ?install.packages")
  }

  ## Check input is a 3d array
  if (length(dim(cij)) != 3) {
    stop("Input 'cij' should be a NxNxS 3d array, where N is the number of atoms and S is the number of simulations")
  }
  ## Check input is built of simmetric matrices
  if (dim(cij)[1] != dim(cij)[2]) {
    stop("Input 'cij' should be a NxNxS 3d array, where N is the number of atoms. The matrix uploaded is not NxNxS")
  }
  if (cutoff.sims > dim(cij)[3]) {
    stop("The cutoff.sims number is greater than the number of simulations in the cij matrix")
  }

  ## Filter by cutoff.cij and sum across simulations
  cut.cij.inds <- (abs(cij) < cutoff.cij)
  count <- array(NA, dim = dim(cij))
  count[!cut.cij.inds] = 1
  cij.sum <- apply(count, c(1:2), sum, na.rm = TRUE)

  ## Mask cij values below cutoff and average across simulations
  cij[cut.cij.inds] = NA
  cij.ave <- apply(cij, c(1:2), mean, na.rm = TRUE)

  ## Mask average values if below cutoff.sims
  cut.sims.inds <- (cij.sum < cutoff.sims)
  cij.ave[cut.sims.inds] = 0 ## Could use NA here
  
  class(cij.ave) = c("dccm", "matrix")
  return(cij.ave)
}
