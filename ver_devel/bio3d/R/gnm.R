#' Gaussian Network Model
#'
#' Perform Gaussian network model (GNM) based normal mode analysis (NMA) for 
#' a protein structure.
#'
#' @details This function builds a Gaussian network model (an isotropic elastic network 
#'    model) for C-alpha atoms and performs subsequent normal mode analysis (NMA). 
#'    The model employs a distance cutoff for the network construction: Atom pairs with 
#'    distance falling  within the cutoff have a harmonic interaction with a uniform force constant; 
#'    Otherwise atoms have no interaction. Output contains N-1 (N, the number of residues) 
#'    non-trivial modes (i.e. the degree of freedom is N-1), which can then be used to 
#'    calculate atomic fluctuations and covariance.  
#'
#' @param x an object of class \code{pdb} as obtained from function \code{\link{read.pdb}}. 
#' @param inds atom and xyz coordinate indices obtained from \code{\link{atom.select}} that 
#'    selects the elements of \code{pdb} upon which the calculation should be based. 
#'    If not provided the function will attempt to select all calpha atoms automatically.
#' @param temp numerical, temperature for which the amplitudes for scaling the atomic 
#'    displacement vectors are calculated. Set \sQuote{temp=NULL} to avoid scaling. 
#' @param keep numerical, final number of modes to be stored. Note that all subsequent analyses 
#'    are limited to this subset of modes. This option is useful for very large structures and 
#'    cases where memory may be limited.
#' @param outmodes atom indices as obtained from \code{\link{atom.select}} specifying the atoms 
#'    to include in the resulting mode object. 
#' @param gamma numerical, global scale of the force constant.
#' @param cutoff numerical, distance cutoff for pair-wise interactions.
#' @param check.connect logical, if TRUE check chain connectivity.
#'
#' @return Returns an object of class \sQuote{gnm} with the following components:
#'    \item{force.constants}{ numeric vector containing the force constants corresponding 
#'       to each mode. }
#'    \item{fluctuations}{ numeric vector of atomic fluctuations. }
#'    \item{U}{ numeric matrix with columns containing the raw eigenvectors. }
#'    \item{L}{ numeric vector containing the raw eigenvalues. }
#'    \item{xyz}{ numeric matrix of class \code{xyz} containing the Cartesian coordinates 
#'       in which the calculation was performed. }
#'    \item{temp}{ numerical, temperature for which the amplitudes for scaling the atomic 
#'       displacement vectors are calculated. }
#'    \item{triv.modes}{ number of trivial modes. }
#'    \item{natoms}{ number of C-alpha atoms. }
#'    \item{call}{ the matched call. }
#'
#' @seealso 
#'      \code{\link{gnm.pdbs}}
#'
#' @author Xin-Qiu Yao & Lars Skjaerven
#' 
#' @references 
#'    Bahar, I. et al. (1997) \emph{Folding Des.} \bold{2}, 173. 
#' 
#' @examples
#'    ## Fetch stucture
#'    pdb <- read.pdb( system.file("examples/1hel.pdb", package="bio3d") )
#'    
#'    ## Calculate normal modes
#'    modes <- gnm(pdb)
#'    
#'    ## Print modes
#'    print(modes)
#'    
#'    ## Plot modes
#'    plot(modes)
#'
gnm <- function(x, ...) {
  UseMethod("gnm")
}

#' @rdname gnm
gnm.pdb <- function(x, inds=NULL, temp=300.0, keep=NULL, outmodes=NULL, 
                      gamma=1.0, cutoff=8.0, check.connect=TRUE, ...) {
 
  pdb <- x
 
  ## Log the call
  cl <- match.call()

  if(!is.pdb(pdb))
    stop("please provide a 'pdb' object as obtained from 'read.pdb()'")

  if(!is.null(outmodes) & !is.select(outmodes))
    stop("provide 'outmodes' as obtained from function atom.select()")
  
  ## Prepare PDB
  ## Take only first frame of multi model PDB files
  if(nrow(pdb$xyz)>1) {
    warning("multimodel PDB file detected - using only first frame")
    pdb$xyz=pdb$xyz[1,, drop=FALSE]
  }

  ## Trim to only CA atoms
  if(is.null(inds)) {
    ca.inds <- atom.select(pdb, "calpha", verbose=FALSE)
    pdb.in <- trim.pdb(pdb, ca.inds)
  }

  ## or to user selection
  else {
    pdb.in <- trim.pdb(pdb, inds)
    if(!all(pdb.in$atom$elety=="CA"))
      stop("non-CA atoms detected")
  }

  ## Indices for effective hessian
  if(is.select(outmodes)) {
    ## re-select since outmodes indices are based on input PDB
    inc.inds <- .match.sel(pdb, pdb.in, outmodes)
    pdb.out <- trim.pdb(pdb.in, inc.inds)
  }
  else {
    pdb.out <- pdb.in
    inc.inds <- atom.select(pdb.in, "all", verbose=FALSE)
  }

  ## fetch number of atoms and sequence
  natoms.in  <- ncol(pdb.in$xyz)/3
  natoms.out <- ncol(pdb.out$xyz)/3
  
  if (natoms.in<2)
    stop("gnm: insufficient number of atoms")

  ## check structure connectivity
  if(check.connect) {
    conn <- inspect.connectivity(pdb.in$xyz)
    if(!conn) {
      warning("Possible multi-chain structure or missing in-structure residue(s) present\n", 
              "  Fluctuations at neighboring positions may be affected.")
    }
  }
 
  ## build Kirchhoff matrix
  K <- cmap(pdb.in, scut=1, dcut=cutoff, mask.lower=FALSE)
  K[!is.na(K) & K > 0] <- -1
  diag(K) <- - apply(K, 1, sum, na.rm=TRUE)
  
  if(length(inc.inds$atom) < nrow(K)) {
    # calculate effective Kirchhoff
    ptm <- proc.time()
#    cat(" Extracting effective Kirchhoff..")
    inc.inds <- inc.inds$atom
    kaa    <- K[inc.inds, inc.inds]
    kqq.inv <- chol2inv(chol(K[-inc.inds, -inc.inds]))
    kaq     <- K[inc.inds, -inc.inds]
    kqa     <- t(kaq)
    K <- kaa - crossprod(crossprod(kqq.inv, kqa), kqa)
    t <- proc.time() - ptm
#    cat("\tDone in", t[[3]], "seconds.\n")
  }

  ## diagonalize - get eigenvectors
  ei <- eigen(K, symmetric=TRUE)
  ei$values <- ei$values[length(ei$values):1]
  ei$vectors <- ei$vectors[, length(ei$values):1]

  if(any(ei$values[-1] < 0)) {
    warning("Negative eigenvalue(s) detected! \
            This can be an indication of an unphysical input structure.")
  }

  ## keep only a subset of modes - including trivial modes
   if(!is.null(keep)) {
    if(keep > ncol(ei$vectors))
      keep <- ncol(ei$vectors)
    keep.inds <- seq(1, keep)
    ei$vectors <- ei$vectors[, keep.inds]
    ei$values <- ei$values[keep.inds]
  }
  
  ## fluctuations
  kb <- 0.00831447086363271
  f <- t(ei$vectors[, -1]**2) / ei$values[-1]
  f <- colSums(f) * 3 * kb * temp / gamma

  ## make a GNM object
  gnm <- list(force.constants = ei$values, 
              fluctuations = f,
              U = ei$vectors, L = ei$values,
              xyz = pdb.out$xyz,
              temp = temp, 
              triv.modes=1,
              natoms = natoms.out,
              call = cl)
  class(gnm) <- c('EnergeticModes', 'gnm', 'nma')
  return(gnm)
}

