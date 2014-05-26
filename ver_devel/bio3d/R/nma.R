
## all-atom NMA
"aanma" <- function(pdb, inds=NULL, ff='aaenm', pfc.fun=NULL, mass=TRUE,
                    temp=300.0, keep=NULL, hessian=NULL, outmodes='calpha', ... ) {
  outmodes.allowed <- c("calpha", "noh")
  if(!outmodes%in%outmodes.allowed)
    stop("outmodes must be 'calpha' or 'noh'")
  
  ## Log the call
  cl <- match.call()
  
  if(!is.pdb(pdb))
    stop("nma: 'pdb' must be of type 'pdb'")
  
  ## Initialize
  init <- .nma.init(ff=ff, pfc.fun=pfc.fun, ...)
  
  ## Trim PDB to match user selection
  if ( !is.null(inds) ) {
    pdb <- trim.pdb(pdb, inds)
  }
  
  aa.inds <- atom.select(pdb, "all", verbose=FALSE)
  ha.inds <- atom.select(pdb, "noh", verbose=FALSE)
  ca.inds <- atom.select(pdb, "calpha", verbose=FALSE)
  ha.pdb <- trim.pdb(pdb, ha.inds)
  ca.pdb <- trim.pdb(pdb, ca.inds)
  
  sequ    <- pdbseq(ca.pdb)
  ca.atoms  <- length(ca.inds$atom)
  ha.atoms  <- length(ha.inds$atom)
  aa.atoms  <- length(aa.inds$atom)
  
  if (ca.atoms<3)
    stop("nma: insufficient number of CA atoms in structure")
  
  ## Use aa2mass to fetch residue mass
  if (mass) {
    aa.masses <- atom2mass(pdb)
    ha.masses <- atom2mass(ha.pdb)
    ca.masses <- atom2mass(ha.pdb, grpby=ha.pdb$atom[,"resno"])
  }
  
  ## Residue mass is provided by user
  else if (!is.null(init$bh.args$aa.mass)) {
    warning("provided masses not in effect at the moment")
  }
  
  ## No mass-weighting
  else {
    ha.masses <- NULL; ca.masses <- NULL;
    init$bh.args <- init$bh.args[ !('aa.mass' %in% names(bh.args)) ]
  }

  ## use 'in' for build.hessian, 'out' for nma.finalize
  masses.in <- ha.masses
  xyz.in <- ha.pdb$xyz
  inc.inds <- NULL
  
  if(!is.null(hessian)) {
    xyz.in <- NULL
    masses.in <- NULL

    if(outmodes=="calpha")
      inc.inds <- atom.select(pdb, "calpha", verbose=FALSE)
    else
      inc.inds <- atom.select(pdb, "noh", verbose=FALSE)
  }

  if(outmodes=="calpha") {
    inc.inds <- atom.select(ha.pdb, "calpha", verbose=FALSE)
    masses.out <- ca.masses
    natoms <- length(ca.masses)
    xyz.out <- ca.pdb$xyz
  }
  else {
    masses.out <- ha.masses
    natoms <- length(ha.masses)
    xyz.out <- ha.pdb$xyz
  }
  
  ## NMA hessian
  ## sequ is here of length ca.pdb, which is different from natoms !!!
  hessian <- .nma.hess(xyz.in, init=init, sequ=sequ, masses=masses.in,
                       hessian=hessian, inc.inds=inc.inds)
  
  ## diagaonalize - get eigenvectors
  ei <- .nma.diag(hessian)

  ## make a NMA object
  nmaobject <- .nma.finalize(ei, xyz=xyz.out, temp=temp, masses=masses.out,
                             natoms=natoms, keep=keep, call=cl)
  return(nmaobject)
}

## calpha NMA
"nma" <- function(pdb, inds=NULL, ff='calpha', pfc.fun=NULL, mass=TRUE,
                  temp=300.0, keep=NULL, hessian=NULL,  ... ) {
  ## Log the call
  cl <- match.call()

  if(!is.pdb(pdb))
    stop("nma: 'pdb' must be of type 'pdb'")
  
  ## Initialize
  init <- .nma.init(ff=ff, pfc.fun=pfc.fun, ...)
  
  ## Trim PDB to match user selection
  if ( !is.null(inds) ) {
    pdb <- trim.pdb(pdb, inds)
  }
  
  ## Trim to only CA atoms
  ca.inds <- atom.select(pdb, "calpha", verbose=FALSE)
  pdb <- trim.pdb(pdb, ca.inds)
  sequ    <- pdbseq(pdb)
  natoms  <- length(ca.inds$atom)
  if (natoms<3)
    stop("nma: insufficient number of CA atoms in structure")
    
  ## Use aa2mass to fetch residue mass
  if (mass && is.null(init$bh.args$aa.mass) ) {
    masses <- do.call('aa2mass', c(list(pdb=sequ, inds=NULL), init$am.args))
  }
  
  ## Residue mass is provided by user
  else if (!is.null(init$bh.args$aa.mass)) {
    masses <- init$bh.args$aa.mass
    init$bh.args <- init$bh.args[ !('aa.mass' %in% names(init$bh.args)) ]
    if(!mass) {
      warning("incompatible arguments: forcing mass weighting")
      mass <- TRUE
    }
  }
  
  ## No mass-weighting
  else {
    masses <- NULL
    init$bh.args <- init$bh.args[ !('aa.mass' %in% names(bh.args)) ]
  }
  
  ## NMA hessian
  hessian <- .nma.hess(pdb$xyz, init=init, sequ=sequ, masses=masses,
                       hessian=hessian, inc.inds=NULL)

  ## diagaonalize - get eigenvectors
  ei <- .nma.diag(hessian)

  ## make a NMA object
  nmaobject <- .nma.finalize(ei, xyz=pdb$xyz, temp=temp, masses=masses,
                             natoms=natoms, keep=keep, call=cl)
  return(nmaobject)
}


".nma.init" <- function(ff=NULL, pfc.fun=NULL, ...) {
  
    ## Arguments to functions build.hessian and aa2mass
    bh.names <- names(formals( build.hessian ))
    am.names <- names(formals( aa2mass ))
    
    dots <- list(...)
    bh.args <- dots[names(dots) %in% bh.names]
    am.args <- dots[names(dots) %in% am.names]   

    ## Define force field
    if (is.null(pfc.fun)) {
      ff <- load.enmff(ff)
    }
    else {
      ## Use customized force field
      if(!is.function(pfc.fun))
        stop("'pfc.fun' must be a function")
      bh.args <- bh.args[ !('pfc.fun' %in% names(bh.args)) ]
      ff <- pfc.fun
    }

    ## Check for optional arguments to pfc.fun
    ff.names <- names(formals( ff ))
    ff.args  <- dots[names(dots) %in% ff.names]

    ## Redirect them to build.hessian
    bh.args  <- c(bh.args, ff.args)

    ## Arguments without destination
    all.names <- unique(c(bh.names, am.names, ff.names))
    if(!all(names(dots) %in% all.names)) {
      oops <- names(dots)[!(names(dots) %in% all.names)]
      stop(paste("argument mismatch:", oops))
    }

    if(length(bh.args)==0)
      bh.args=NULL
    if(length(am.args)==0)
      am.args=NULL
    if(length(ff.args)==0)
      ff.args=NULL
        
    out <- list(pfcfun=ff, bh.args=bh.args, ff.args=ff.args, am.args=am.args)
    return(out)
  }



## effective hessian
".nma.trim.hessian" <- function(hess, inc.inds=NULL) {
  if(!is.matrix(hess))
    stop("hess must be a matrix")
  if(is.null(inc.inds))
    stop("indices must be provided")
  
  kaa     <- hess[inc.inds, inc.inds]
  kqq.inv <- solve(hess[-inc.inds, -inc.inds])
  kaq     <- hess[inc.inds, -inc.inds]
  kqa     <- t(kaq)
  k <- kaa - ((kaq %*% kqq.inv) %*% kqa)
  return(k)
}

".nma.hess" <- function(xyz, init, sequ, masses, hessian, inc.inds=NULL) {
  
  ## Build the Hessian Matrix
  if(is.null(hessian)) {
    cat(" Building Hessian...")
    ptm <- proc.time()
    H <- do.call('build.hessian', list(xyz=xyz, pfc.fun=init$pfcfun, sequ=sequ, aa.mass=masses))
    ##H <- build.hessian(xyz, pfc.fun=pfcfun, aa.mass=masses)
    t <- proc.time() - ptm
    cat("\t\tDone in", t[[3]], "seconds.\n")
  }
  else {
    H <- hessian
  }
  
  ## Effective Hessian
  if(!is.null(inc.inds)) {
    cat(" Extracting effective Hessian...")
    ptm <- proc.time()
    H <- .nma.trim.hessian(H, inc.inds=inc.inds$xyz)
    t <- proc.time() - ptm
    cat("\t\tDone in", t[[3]], "seconds.\n")
  }
  return(H)
}

  
## build a NMA object
".nma.diag" <- function(H) {
  
  ## Diagonalize matrix
  cat(" Diagonalizing Hessian...")
  ptm <- proc.time()
  ei <- eigen(H, symmetric=TRUE)
  t <- proc.time() - ptm
  cat("\tDone in", t[[3]], "seconds.\n")
  return(ei)
}

".nma.finalize" <- function(ei, xyz, temp, masses, natoms, keep, call) {
  if(length(masses)>0)
    mass <- TRUE
  else
    mass <- FALSE
  
  if(!is.null(keep)) {
    if(keep>ncol(ei$vectors))
      keep <- ncol(ei$vectors)
    keep <- keep-1
    keep.inds <- seq(ncol(ei$vectors)-keep, ncol(ei$vectors))
    ei$vectors <- ei$vectors[,keep.inds]
    ei$values <- ei$values[keep.inds]
  }
  
  ## Raw eigenvalues
  ei$values <- round(ei$values, 6)
  ##triv.modes <- which(ei$values<=0) ## indicies !!

  ## hard code 6 trivial modes instead
  triv.modes <- seq(length(ei$values), length(ei$values)-5)

    ## Frequencies are given by
    if (mass)  {
      pi <- 3.14159265359
      freq <- sqrt(abs(ei$values)) / (2 * pi)
      force.constants <- NULL
    } else {
      freq <- NULL
      force.constants <- ei$values
    }

    ## Raw unmodified eigenvectors:
    ## ei$vectors

    ## V holds the eigenvectors converted to unweighted Cartesian coords:
    V <- ei$vectors

    ## Change to non-mass-weighted eigenvectors
    if(mass) {
      wts.sqrt <- sqrt(masses)
      tri.inds <- rep(1:natoms, each=3)
      V <- apply(V, 2, '*', 1 / wts.sqrt[tri.inds])
    }

    ## Temperature scaling
    kb <- 0.00831447086363271
    if ( !is.null(temp) ) {
      if (!is.null(freq)) {
        amplitudes <- sqrt(2* temp * kb) / (2* pi * freq[ -triv.modes ])
        amplitudes <- c(amplitudes, rep(1,length(triv.modes)))
      }
      else if(!is.null(force.constants)) {
        amplitudes <- sqrt((2* temp * kb) /
                           force.constants[ -triv.modes ])
        amplitudes <- c(amplitudes, rep(1,length(triv.modes)))
      }
    } else {
      amplitudes <- rep(1, times=3*natoms)
    }

    ## Trivial modes first (reverse matrix!)
    ei$vectors <- ei$vectors[, seq(ncol(ei$vectors),1)]
    V <- V[, seq(ncol(ei$vectors),1)]
    ei$values <- rev( ei$values )
    freq <- rev(freq)
    force.constants <- rev(force.constants)
    amplitudes <- rev(amplitudes)

    ## Temperature scaling of eigenvectors
    for ( i in (length(triv.modes)+1):ncol(V) ) {
      V[,i] <- (V[,i] * amplitudes[i])
    }

    ## Check if first modes are zero-modes
    if(ei$values[1]<0) {
      warning("Negative eigenvalue(s) detected! \
              This can be an indication of an unphysical input structure.")
    }

    ## Output to class "nma"
    nma <- list(modes=V,
                frequencies=NULL,
                force.constants=NULL,
                fluctuations=NULL,
                U=ei$vectors, L=ei$values,

                xyz=xyz,
                mass=masses,
                temp=temp,
                triv.modes=length(triv.modes),
                natoms=natoms,
                call=call)

    if(mass) {
      class(nma) <- c("VibrationalModes", "nma")
      nma$frequencies <- freq
    }
    else {
      class(nma) <- c("EnergeticModes", "nma")
      nma$force.constants <- force.constants
    }

    ## Calculate mode fluctuations
    nma$fluctuations <- fluct.nma(nma, mode.inds=NULL)

    ## Notes:
    ## U are the raw unmodified eigenvectors
    ## These mode vectors are in mass-weighted coordinates and not
    ## scaled by the thermal amplitudes, so they are orthonormal.

    ## V holds the eigenvectors converted to unweighted Cartesian
    ## coordinates.  Unless you set temp=NULL, the modes are
    ## also scaled by the thermal fluctuation amplitudes.

    return(nma)
  }
