".nma.args" <- function(pfc.fun=NULL, ...) {

    ## Arguments to functions build.hessian and aa2mass
    bh.names <- names(formals( build.hessian ))
    am.names <- names(formals( aa2mass ))
    #rtb.names <- names(formals( rtb ))

    dots <- list(...)
    bh.args <- dots[names(dots) %in% bh.names]
    bh.args <- bh.args[ !('pfc.fun' %in% names(bh.args)) ]
    am.args <- dots[names(dots) %in% am.names]
    #rtb.args <- dots[names(dots) %in% rtb.names]

    ## Check for optional arguments to pfc.fun
    ff.names <- names(formals( pfc.fun ))
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
    #if(length(rtb.args)==0)
    #    rtb.args=NULL

    out <- list(bh.args=bh.args, am.args=am.args) #, rtb.args=rtb.args)
    return(out)
  }

## extract effective hessian for RTB approach
.nma.trim.hessian.rtb <- function(hessian, inc.inds=NULL, pdb=NULL, nmer=1) {

    if(!is.matrix(hessian))
      stop("hessian must be a matrix")
    if(is.null(inc.inds))
      return(hessian)

    if(nrow(hessian) == length(inc.inds$xyz))
      return(hessian)

    cat(" Extracting effective Hessian with RTB..")
    ptm <- proc.time()

    exc.inds <- combine.select(atom.select(pdb),
                               inc.inds, operator='-', verbose=FALSE)
    pdb <- trim(pdb, exc.inds)

    inc.inds <- inc.inds$xyz
    kaa <- hessian[inc.inds, inc.inds]
    ei <- rtb(hessian[-inc.inds, -inc.inds], pdb=pdb, mass=FALSE, nmer=nmer,
              verbose=FALSE)

    if(any(ei$values<=0)) {
      warning("Not all eigenvalues are positive!")
    }
    inds <- which(ei$values > 0)

    kaq     <- hessian[inc.inds, -inc.inds]
    kqa     <- t(kaq)
    ukqa <- crossprod(ei$vectors[, inds], kqa)
    k <- kaa - crossprod(ukqa * (1/ei$values[inds]), ukqa)

    t <- proc.time() - ptm
    cat("\tDone in", t[[3]], "seconds.\n")

    return(k)
}


## extract effective hessian
".nma.trim.hessian" <- function(hessian, inc.inds=NULL) {

  if(!is.matrix(hessian))
    stop("hessian must be a matrix")
  if(is.null(inc.inds))
    return(hessian)

  if(nrow(hessian) == length(inc.inds$xyz))
    return(hessian)

  ptm <- proc.time()
  cat(" Extracting effective Hessian..")
  inc.inds <- inc.inds$xyz

  kaa    <- hessian[inc.inds, inc.inds]
  ##kqq.inv <- solve(hessian[-inc.inds, -inc.inds])
  kqq.inv <- chol2inv(chol(hessian[-inc.inds, -inc.inds]))
  kaq     <- hessian[inc.inds, -inc.inds]
  kqa     <- t(kaq)
  ##k <- kaa - ((kaq %*% kqq.inv) %*% kqa)
  k <- kaa - crossprod(crossprod(kqq.inv, kqa), kqa)

  t <- proc.time() - ptm
  cat("\tDone in", t[[3]], "seconds.\n")
  return(k)
}

## mass-weight hessian
".nma.mwhessian" <- function(hessian, masses=NULL) {
  if(!is.matrix(hessian))
    stop("hessian must be a matrix")
  if(is.null(masses))
    stop("masses must be provided")

  #cat(" Mass weighting Hessian...")
  #ptm <- proc.time()

  dims <- dim(hessian)
  natoms <- dims[1] / 3

  if(length(masses)!=natoms)
    stop("dimension mismatch")

  masses <- sqrt(masses)
  inds <- rep(1:natoms, each=3)
  col.inds <- seq(1, ncol(hessian), by=3)

  for ( i in 1:natoms ) {
    m <- col.inds[i]
    hessian[,m:(m+2)] <- hessian[,m:(m+2)] * (1/masses[i])
    hessian[,m:(m+2)] <- hessian[,m:(m+2)] * (1/masses[inds])
  }

  #t <- proc.time() - ptm
  #cat("\tDone in", t[[3]], "seconds.\n")

  return(hessian)
}


## wrapper for generating the hessian matrix
".nma.hess" <- function(xyz, pfc.fun, args=NULL,
                        hessian=NULL, pdb=NULL) {

  natoms <- ncol(as.xyz(xyz))/3
  if(nrow(xyz)>1)
    xyz=xyz[1,,drop=FALSE]

  ## Build the Hessian Matrix
  if(is.null(hessian)) {
    cat(" Building Hessian...")
    ptm <- proc.time()
    H <- do.call('build.hessian',
                 c(list(xyz=xyz, pfc.fun=pfc.fun, pdb=pdb), args$bh.args))
    t <- proc.time() - ptm
    cat("\t\tDone in", t[[3]], "seconds.\n")
  }
  else {
    H <- hessian
  }
  return(H)
}


## diagonalize hessian
".nma.diag" <- function(H, symmetric=TRUE) {

  ## Diagonalize matrix
  cat(" Diagonalizing Hessian...")
  ptm <- proc.time()
  ei <- eigen(H, symmetric=symmetric)
  t <- proc.time() - ptm
  cat("\tDone in", t[[3]], "seconds.\n")
  return(ei)
}

## build a NMA object
".nma.finalize" <- function(ei, xyz, temp, masses, natoms, keep, call) {
  if(length(masses)>0)
    mass <- TRUE
  else
    mass <- FALSE

  xyz=as.xyz(xyz)
  dims <- dim(ei$vectors)
  dimchecks <- c(ncol(xyz)/3==natoms,
                 ifelse(mass, length(masses)==natoms, TRUE),
                 dims[1]/3==natoms)
                 #dims[2]/3==natoms)

  if(!all(dimchecks))
    stop(paste("dimension mismatch when generating nma object\n",
               paste(dimchecks, collapse=", ")))

  ## Raw eigenvalues
  ei$values <- round(ei$values, 6)

  ## Trivial modes first - sort on abs(ei$values)
  sort.inds  <- order(abs(ei$values))
  ei$values  <- ei$values[sort.inds]
  ei$vectors <- ei$vectors[, sort.inds]

  ## hard code 6 trivial modes
  triv.modes <- seq(1, 6)

  ## keep only a subset of modes - including trivial modes
   if(!is.null(keep)) {
    if(keep>ncol(ei$vectors))
      keep <- ncol(ei$vectors)
    keep.inds <- seq(1, keep)
    ei$vectors <- ei$vectors[,keep.inds]
    ei$values <- ei$values[keep.inds]
  }

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
        amplitudes <- c(rep(1,length(triv.modes)), amplitudes)
      }
      else if(!is.null(force.constants)) {
        amplitudes <- sqrt((2* temp * kb) /
                           force.constants[ -triv.modes ])
        amplitudes <- c(rep(1,length(triv.modes)), amplitudes)
      }
    } else {
      amplitudes <- rep(1, times=3*natoms)
    }

    ## Temperature scaling of eigenvectors
    for ( i in (length(triv.modes)+1):ncol(V) ) {
      V[,i] <- (V[,i] * amplitudes[i])
    }

    ## Check if first modes are zero-modes
    if(any(ei$values<0)) {
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


.inds2ids <- function(pdb, inds=NULL) {
    if(!is.null(inds)) {
        paste(pdb$atom$chain[inds$atom],
              pdb$atom$resno[inds$atom],
              pdb$atom$elety[inds$atom],
              pdb$atom$insert[inds$atom], sep="_")
    }
    else {
        paste(pdb$atom$chain,
              pdb$atom$resno,
              pdb$atom$elety,
              pdb$atom$insert, sep="_")
    }
}

".match.sel" <- function(a, b, inds) {
  ## a= original pdb
  ## b= trimmed pdb
  ## inds= indices of pdb 'a' to keep
  ## find corresponding atoms in b

  ids.a <- .inds2ids(a, inds)
  ids.b <- .inds2ids(b)

  inds <- which(ids.b %in% ids.a)
  return(as.select(inds))
}


".nma.reduce.pdb" <- function(pdb) {

    pdb$atom$resid <- aa123(aa321(pdb$atom$resid))

    ## Selected side chain atoms
    ## Using whatever is furthest away from CA, including N and O atoms.
    aa.elety <- list(
        ARG = c("NH1", "NH2"),
        HIS = c("CE1", "CG"),
        LYS = c("NZ"),
        ASP = c("OD1", "OD2"),
        GLU = c("OE1", "OE2"),
        SER = c("OG"),
        THR = c("CG2", "OG1"), ## "OG1"
        ASN = c("OD1", "ND2"),
        GLN = c("NE2", "OE1"),
        CYS = c("SG"),
        ##GLY = c(),
        PRO = c("CG"),
        ALA = c("CB"),
        VAL = c("CG1", "CG2"),
        ILE = c("CD1", "CG2"),
        LEU = c("CD1", "CD2"),
        MET = c("CE"),
        PHE = c("CZ"),
        TYR = c("OH"),
        TRP = c("CH2", "NE1")
        )

    pdb <- trim(pdb, "noh")
    sele <- atom.select(pdb, elety=c("N", "CA", "C"))

    for(i in 1:length(aa.elety)) {
        resid = names(aa.elety)[i]
        elety = aa.elety[[i]]
        sele2 <- atom.select(pdb, resid = resid, elety = elety)

        if(length(sele2$atom)>0) {
            sele <- combine.select(sele, sele2, operator="OR", verbose=FALSE)
        }
    }

    #return(sele)
    return(trim(pdb, sele))

}
