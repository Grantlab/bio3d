

## all-atom NMA
"aanma" <- function(pdb, ff='aaenm', pfc.fun=NULL, mass=TRUE,
                    temp=300.0, keep=NULL, hessian=NULL, outmodes='calpha', ... ) {
  
  ## Log the call
  cl <- match.call()

  ## Indices for effective hessian
  if(!outmodes %in% c("calpha", "noh") && !is.select(outmodes) && is.null(hessian))
    stop("outmodes must be 'calpha', 'noh', or an atom selection by 'atom.select()'")
      
  if(!is.pdb(pdb))
    stop("nma: 'pdb' must be of type 'pdb'")
  
  ## Initialize
  init <- .nma.init(ff=ff, pfc.fun=pfc.fun, ...)

  if(!is.null(hessian)) {
    pdb.in <- pdb
    dims <- dim(hessian)
    if(dims[1]!=dims[2] | dims[1]!=length(pdb.in$xyz))
      stop("dimension mismatch")
  }
  else {
    tmp.inds <- atom.select(pdb, "noh", verbose=FALSE)
    pdb.in <- trim.pdb(pdb, tmp.inds)
  }
  
  ## Indices
  if(is.select(outmodes)) {
    ## since pdb.in is 'noh' (from trim.pdb above), we need to re-select :P
    tmp.inds <- outmodes
    inc.inds <- .match.sel(pdb, pdb.in, tmp.inds)
    pdb.out <- trim.pdb(pdb.in, inc.inds)
  }
  else if(outmodes=="calpha") {
    inc.inds <- atom.select(pdb.in, "calpha", verbose=FALSE)
    pdb.out <- trim.pdb(pdb.in, inc.inds)
  }
  else {
    inc.inds <- atom.select(pdb.in, "noh", verbose=FALSE)
    pdb.out <- trim.pdb(pdb.in, inc.inds)
  }
  
  sequ    <- pdbseq(pdb.in)
  natoms.in <- length(pdb.in$xyz)/3
  natoms.out <- length(pdb.out$xyz)/3
  
  if (natoms.in<3 | natoms.out<3)
    stop("nma: insufficient number of atoms in selection")
  
  ## Use aa2mass to fetch residue mass
  if (mass) {
    masses.in <-  atom2mass(pdb.in)
    ##masses.out <- atom2mass(pdb.out, grpby=pdb.out$atom[,"resno"])
    masses.out <- aa2mass(pdb.out)
  }
  
  ## Residue mass is provided by user
  #else if (!is.null(init$bh.args$aa.mass)) {
  #  warning("provided masses not in effect at the moment")
  #}
  
  ## No mass-weighting
  else {
    masses.in <- NULL; masses.out <- NULL;
    init$bh.args <- init$bh.args[ !('aa.mass' %in% names(bh.args)) ]
  }
  
  ## NMA hessian
  ## sequ is here of length ca.pdb, which is different from natoms !!!
  hessian <- .nma.hess(pdb.in$xyz, init=init, sequ=sequ, masses=masses.in,
                       hessian=hessian, inc.inds=inc.inds)
  
  ## diagaonalize - get eigenvectors
  ei <- .nma.diag(hessian)

  ## make a NMA object
  nmaobject <- .nma.finalize(ei, xyz=pdb.out$xyz, temp=temp, masses=masses.out,
                             natoms=natoms.out, keep=keep, call=cl)
  return(nmaobject)
}

## calpha NMA
"nma" <- function(pdb, ff='calpha', pfc.fun=NULL, mass=TRUE,
                  temp=300.0, keep=NULL, hessian=NULL, outmodes=NULL, ... ) {
  
  ## Log the call
  cl <- match.call()

  if(!is.pdb(pdb))
    stop("nma: 'pdb' must be of type 'pdb'")
  
  ## Initialize
  init <- .nma.init(ff=ff, pfc.fun=pfc.fun, ...)
  
  ## Trim to only CA atoms
  ca.inds <- atom.select(pdb, "calpha", verbose=FALSE)
  pdb.in <- trim.pdb(pdb, ca.inds)
  
  if(is.select(outmodes)) {
    ## since pdb.in is 'noh' (from trim.pdb above), we need to re-select :P
    tmp.inds <- outmodes
    inc.inds <- .match.sel(pdb, pdb.in, tmp.inds)
    pdb.out <- trim.pdb(pdb.in, inc.inds)
  }
  else {
    pdb.out <- pdb.in
    inc.inds <- atom.select(pdb.in, "all")
  }
  
  sequ    <- pdbseq(pdb.in)
  natoms.in <- length(pdb.in$xyz)/3
  natoms.out <- length(pdb.out$xyz)/3

  if (natoms.in<3)
    stop("nma: insufficient number of CA atoms in structure")
  
  ## Use aa2mass to fetch residue mass
  ##if (mass && is.null(init$bh.args$aa.mass) ) {
  if (mass) {
    masses.in <- do.call('aa2mass', c(list(pdb=sequ, inds=NULL), init$am.args))
    masses.out <- masses.in[ inc.inds$atom ]
  }
  
  ## Residue mass is provided by user
  #else if (!is.null(init$bh.args$aa.mass)) {
  #  masses <- init$bh.args$aa.mass
  #  init$bh.args <- init$bh.args[ !('aa.mass' %in% names(init$bh.args)) ]
  #  if(!mass) {
  #    warning("incompatible arguments: forcing mass weighting")
  #    mass <- TRUE
  #  }
  #}
  
  ## No mass-weighting
  else {
    masses.in <- NULL; masses.out <- NULL;
    init$bh.args <- init$bh.args[ !('aa.mass' %in% names(bh.args)) ]
  }
  
  ## NMA hessian
  hessian <- .nma.hess(pdb.in$xyz, init=init, sequ=sequ, masses=masses.in,
                       hessian=hessian, inc.inds=inc.inds)

  ## diagaonalize - get eigenvectors
  ei <- .nma.diag(hessian)

  ## make a NMA object
  nmaobject <- .nma.finalize(ei, xyz=pdb.out$xyz, temp=temp, masses=masses.out,
                             natoms=natoms.out, keep=keep, call=cl)
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
    t <- proc.time() - ptm
    cat("\t\tDone in", t[[3]], "seconds.\n")
  }
  else {
    H <- hessian
  }
  
  ## Effective Hessian
  if(!is.null(inc.inds)) {
    if(length(xyz)>length(inc.inds$xyz)) {
      cat(" Extracting effective Hessian..")
      ptm <- proc.time()
      H <- .nma.trim.hessian(H, inc.inds=inc.inds$xyz)
      t <- proc.time() - ptm
      cat("\tDone in", t[[3]], "seconds.\n")
    }
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

".match.sel" <- function(a, b, inds) {
  ## a= original pdb
  ## b= trimmed pdb
  ## inds= indices of pdb 'a' to keep
  ## find corresponding atoms in b
  
  names.a <- paste(a$atom[inds$atom, "resno"],
                   a$atom[inds$atom, "elety"],
                   a$atom[inds$atom, "eleno"], sep="-")
  
  names.b <- paste(b$atom[, "resno"],
                   b$atom[, "elety"],
                   b$atom[, "eleno"], sep="-")
  
  inds <- which(names.b %in% names.a)
  out <- list(atom=inds, xyz=atom2xyz(inds))
  class(out) <- "select"
  return(out)
}
