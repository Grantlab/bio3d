#' All Atom Normal Mode Analysis 
#' 
#' Perform all-atom elastic network model normal modes calculation of a protein 
#' structure. 
#' 
#' @details This function builds an elastic network model (ENM) based on all  
#' heavy atoms of input \code{pdb}, and performs subsequent normal mode 
#' analysis (NMA) in various manners. By default, the \sQuote{aaenm2} force 
#' field (defining of the spring constants between atoms) is used, which was
#' obtained by fitting to a local energy minimum of a crambin model 
#' derived from the AMBER99SB force field. It employs a pair force constant
#' function which falls as r^-6, and specific force constants for
#' covalent and intra-residue atom pairs. See also \code{\link{load.enmff}} 
#' for other force field options. 
#' 
#' The \code{outmodes} argument controls the type of output modes. There are 
#' two standard types of output modes: \sQuote{noh} and \sQuote{calpha}.
#' \code{outmodes='noh'} invokes regular all-atom based ENM-NMA. When
#' \code{outmodes='calpha'}, an effective Hessian with respect to all C-alpha
#' atoms will be first calculated using the same formula as in Hinsen et al.  
#' NMA is then performed on this effective C-alpha based Hessian. In addition,
#' users can provide their own atom selection (see \code{\link{atom.select}}) 
#' as the value of \code{outmodes} for customized output modes generation. 
#' 
#' When \code{reduced=TRUE}, only a selection of all heavy atoms is used 
#' to build the ENM. More specifically, three to five atoms per residue
#' constitute the model. Here the N, CA, C atoms represent the protein
#' backbone, and zero to two selected side chain atoms represent the side chain 
#' (selected based on side chain size and the distance to CA). This 
#' coarse-grained ENM has significantly improved computational efficiency and 
#' similar prediction accuracy with respect to the all-atom ENM.
#'  
#' When \code{rtb=TRUE}, rotation-translation block (RTB) based approximate 
#' modes will be calculated. In this method, each residue is assumed to be a 
#' rigid body (or \sQuote{block}) that has only rotational and translational 
#' degrees of freedom. Intra-residue deformation is thus ignored. 
#' (See Durand et al 1994 and Tama et al. 2000 for more details). N residues per
#' block is also supported, where N=1, 2, 3, etc. (See argument \code{nmer}).  
#' The RTB method has significantly improved computational efficiency and 
#' similar prediction accuracy with respect to the all-atom ENM.
#' 
#' By default the function will diagonalize the mass-weighted Hessian matrix. 
#' The resulting mode vectors are moreover scaled by the thermal fluctuation 
#' amplitudes.
#'
#' @param pdb an object of class \code{pdb} as obtained from function 
#'    \code{\link{read.pdb}}. 
#' @param pfc.fun customized pair force constant (\sQuote{pfc}) function. The 
#'    provided function should take a vector of distances as an argument to 
#'    return a vector of force constants. If NULL, the default function 
#'    \sQuote{aaenm2} will be employed. (See details below).
#' @param mass logical, if TRUE the Hessian will be mass-weighted.
#' @param temp numerical, temperature for which the amplitudes for scaling the 
#'    atomic displacement vectors are calculated. Set \sQuote{temp=NULL} to 
#'    avoid scaling. 
#' @param keep numerical, final number of modes to be stored. Note that all 
#'    subsequent analyses are limited to this subset of modes. This option is 
#'    useful for very large structures and cases where memory may be limited.
#' @param hessian hessian matrix as obtained from \code{\link{build.hessian}}. 
#'    For internal purposes and generally not intended for public use.
#' @param outmodes either a character (\sQuote{calpha} or \sQuote{noh}) or atom 
#'    indices as obtained from \code{\link{atom.select}} specifying the atoms to
#'    include in the resulting mode object. (See details below). 
#' @param rm.wat logical, if TRUE water molecules will be removed before 
#'    calculation. 
#' @param reduced logical, if TRUE the coarse-grained (\sQuote{4-bead}) ENM will
#'    be employed. (See details below).
#' @param rtb logical, if TRUE the rotation-translation block based 
#'    approximate modes will be calculated. (See details below).
#' @param nmer numerical, defines the number of residues per block (used only
#'     when \code{rtb=TRUE}).
#' @param ... additional arguments to \code{\link{build.hessian}} and
#'    \code{\link{aa2mass}}. One useful option here for dealing with 
#'    unconventional residues is \sQuote{mass.custom}, see the 
#'    \code{\link{aa2mass}} function for details.
#' 
#' @return Returns an object of class \sQuote{nma} with the following 
#'    components:
#'    \item{modes}{ numeric matrix with columns containing the normal mode
#'       vectors. Mode vectors are converted to unweighted Cartesian
#'       coordinates  when \code{mass=TRUE}. Note that the 6 first trivial
#'       eigenvectos appear in columns one to six. }
#'    \item{frequencies}{ numeric vector containing the vibrational
#'       frequencies corresponding to each mode (for \code{mass=TRUE}). }
#'    \item{force.constants}{ numeric vector containing the force constants
#'       corresponding to each mode (for \code{mass=FALSE)}). }
#'    \item{fluctuations}{ numeric vector of atomic fluctuations. }
#'    \item{U}{ numeric matrix with columns containing the raw
#'       eigenvectors. Equals to the \code{modes} component when
#'      \code{mass=FALSE} and \code{temp=NULL}. }
#'    \item{L}{ numeric vector containing the raw eigenvalues. }
#'    \item{xyz}{ numeric matrix of class \code{xyz} containing the
#'      Cartesian coordinates in which the calculation was performed. }
#'    \item{mass}{ numeric vector containing the residue masses used for the
#'      mass-weighting. }
#'    \item{temp}{ numerical, temperature for which the amplitudes for
#'      scaling the atomic displacement vectors are calculated. }
#'    \item{triv.modes}{ number of trivial modes. }
#'    \item{natoms}{ number of C-alpha atoms. }
#'    \item{call}{ the matched call. }
#'
#' @seealso 
#'    \code{\link{nma.pdb}} for C-alpha based NMA, \code{\link{aanma.pdbs}} for 
#'    ensemble all-atom NMA, \code{\link{load.enmff}} for available ENM force
#'    fields, and \code{\link{fluct.nma}}, \code{\link{mktrj.nma}}, and 
#'    \code{\link{dccm.nma}} for various post-NMA calculations.
#'
#' @author Lars Skjaerven & Xin-Qiu Yao
#' 
#' @references 
#'    Hinsen, K. et al. (2000) \emph{Chem. Phys.} \bold{261}, 25. 
#'    Durand, P. et al. (1994) \emph{Biopolymers} \bold{34}, 759.
#'    Tama, F. et al. (2000) \emph{Proteins} \bold{41}, 1.
#'    
#' @examples
#' \dontrun{ 
#'    # All-atom NMA takes relatively long time - Don't run by default.
#'  
#'    ## Fetch stucture
#'    pdb <- read.pdb( system.file("examples/1hel.pdb", package="bio3d") )
#'
#'    ## Calculate all-atom normal modes
#'    modes.aa <- aanma(pdb, outmodes='noh')
#'
#'    ## Calculate all-atom normal modes with RTB approximation
#'    modes.aa.rtb <- aanma(pdb, outmodes='noh', rtb=TRUE)
#'
#'    ## Compare the two modes
#'    rmsip(modes.aa, modes.aa.rtb)
#'
#'    ## Calculate C-alpha normal modes.
#'    modes <- aanma(pdb)
#' 
#'    ## Calculate C-alpha normal modes with reduced ENM.
#'    modes.cg <- aanma(pdb, reduced=TRUE)
#'    
#'    ## Calculate C-alpha normal modes with RTB approximation
#'    modes.rtb <- aanma(pdb, rtb=TRUE)
#'
#'    ## Compare modes
#'    rmsip(modes, modes.cg)
#'    rmsip(modes, modes.rtb)
#'
#'    ## Print modes
#'    print(modes)
#'
#'    ## Plot modes
#'    plot(modes)
#'
#'    ## Visualize modes
#'    #m7 <- mktrj.nma(modes, mode=7, file="mode_7.pdb", pdb=pdb)
#' }
"aanma" <- function(pdb, pfc.fun=NULL, mass=TRUE,
                    temp=300.0, keep=NULL, hessian=NULL, outmodes='calpha',
                    rm.wat=TRUE, reduced=FALSE, rtb=FALSE, nmer=1, ... ) {

  ## Log the call
  cl <- match.call()

  if(!outmodes %in% c("calpha", "noh") && !is.select(outmodes))
    stop("outmodes must be 'calpha', 'noh', or an atom selection by 'atom.select()'")

  if(!is.pdb(pdb))
    stop("please provide a 'pdb' object as obtained from 'read.pdb()'")

  ## Define force field
  if (is.null(pfc.fun)) {
      pfc.fun <- load.enmff("aaenm2")
  }
  else {
      ## Use customized force field
      if(!is.function(pfc.fun))
          stop("'pfc.fun' must be a function")
  }

  ## Parse additional arguments
  args <- .nma.args(pfc.fun=pfc.fun, ...)

  ## check and prepare input PDB
  if(!is.null(hessian)) {
    pdb.in <- pdb
    dims <- dim(hessian)
    if(dims[1]!=dims[2] | dims[1]!=length(pdb.in$xyz))
      stop("dimension mismatch")
  }
  else {
    if(rm.wat)
      pdb <- trim.pdb(pdb, "notwater", verbose=FALSE)

    pdb <- trim.pdb(pdb, "noh", verbose=FALSE)
    pdb.in <- pdb

    lig.inds <- atom.select(pdb.in, "ligand")
    if(length(lig.inds$atom)>0) {
      ligs <- paste(unique(pdb.in$atom$resid[ lig.inds$atom ]), sep=", ")
      warning(paste("ligands", ligs, "included in calculation of normal modes"))
    }
  }

  ## reduced aaENM
  if(reduced) {
      pdb.tmp <- .nma.reduce.pdb(pdb.in)

      if(is.select(outmodes)) {
          outmodes <- .match.sel(pdb, pdb.tmp, outmodes)
      }
      pdb.in <- pdb.tmp
      rm(pdb.tmp)
  }

  ## Indices for effective hessian
  ## (selection, calphas, or all (noh) atoms)
  if(is.select(outmodes)) {
    ## since pdb.in is 'noh' (from trim.pdb above), we need to re-select
    inc.inds <- .match.sel(pdb, pdb.in, outmodes)
    pdb.out <- trim.pdb(pdb.in, inc.inds)

    unq.elety <- unique(pdb.out$atom$elety)
    outmodes="noh"
    if(length(unq.elety)==1)
      if(unq.elety=="CA")
        outmodes="calpha"
  }
  else if(outmodes=="calpha") {
    inc.inds <- atom.select(pdb.in, "calpha", verbose=FALSE)
    pdb.out <- trim.pdb(pdb.in, inc.inds)
  }
  else {
    inc.inds <- atom.select(pdb.in, "noh", verbose=FALSE)
    pdb.out <- trim.pdb(pdb.in, inc.inds)
  }

  natoms.in <- nrow(pdb.in$atom)
  natoms.out <- nrow(pdb.out$atom)
  sequ <- pdb.in$atom$resid[ !duplicated(pdb.in$atom$resno) ]

  if (natoms.in<3 | natoms.out<3)
    stop("aanma: insufficient number of atoms")

  ## Masses for weighting the hessian
  ## Note that mass weighting is done in rtb() for rtb=TRUE,
  ## and in .nma.mwhessian() when rtb=FALSE
  if (mass) {
    ## Use atom2mass to fetch atom mass
    if(outmodes=="noh") {
      masses.in <-  atom2mass(pdb.in)
      masses.out <- masses.in[ inc.inds$atom ]
    }

    if(outmodes=="calpha") {
      masses.in <-  atom2mass(pdb.in)
      masses.out <- do.call('aa2mass', c(list(pdb=pdb.out, inds=NULL), args$am.args))
    }
  }
  else {
    ## No mass-weighting
    masses.out <- NULL;
  }

  ## build full hessian
  H <- .nma.hess(pdb.in$xyz, pfc.fun=pfc.fun, args=args,
                       hessian=hessian, pdb=pdb.in)


  ## extract effective hessian and diagonalize
  if(rtb) {
      ## effective hessian
      H <- .nma.trim.hessian.rtb(H, inc.inds=inc.inds, pdb=pdb.in, nmer=nmer)

      ## mass weighting + diagonalize using RTB
      ei <- rtb(H, pdb=pdb.out, mass=masses.out, nmer=nmer)
  }
  else {
      ## effective hessian
      H <- .nma.trim.hessian(H, inc.inds=inc.inds)

      ## mass weight hessian
      if(!is.null(masses.out))
          H <- .nma.mwhessian(H, masses=masses.out)

      ## diagaonalize and obtain eigenvectors
      ei <- .nma.diag(H)
  }

  ## make NMA object
  m <- .nma.finalize(ei, xyz=pdb.out$xyz, temp=temp, masses=masses.out,
                             natoms=natoms.out, keep=keep, call=cl)
  return(m)
}
