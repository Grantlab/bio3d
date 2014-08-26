"nma.pdb" <-
  function(pdb, inds=NULL, ff='calpha', pfc.fun=NULL, mass=TRUE,
           temp=300.0, keep=NULL, hessian=NULL,  ... ) {
    
    if(missing(pdb))
      stop("nma: must supply 'pdb' object, i.e. from 'read.pdb'")
    if(!inherits(pdb, "pdb"))
      stop("nma: 'pdb' must be of type 'pdb'")

    ## Log the call
    cl <- match.call()

    ## Passing arguments to functions build.hessian and aa2mass
    bh.names <- names(formals( build.hessian ))
    am.names <- names(formals( aa2mass ))

    dots <- list(...)
    bh.args <- dots[names(dots) %in% bh.names]
    am.args <- dots[names(dots) %in% am.names]

    ## Trim PDB to match user selection
    if ( !is.null(inds) ) {
      pdb <- trim.pdb(pdb, inds)
    }

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

    ## Only C-alpha ENM NMA is implemented
    ca.inds <- atom.select(pdb, "calpha", verbose=FALSE)
    sequ    <- pdb$atom[ca.inds$atom,"resid"]
    natoms  <- length(ca.inds$atom)
    if (natoms<3)
      stop("nma: insufficient number of CA atoms in structure")
    xyz <- pdb$xyz[ca.inds$xyz]

    ## Use aa2mass to fetch residue mass
    if (mass && is.null(bh.args$aa.mass) ) {
      masses <- do.call('aa2mass', c(list(pdb=sequ, inds=NULL), am.args))
    }

    ## Residue mass is provided by user
    else if (!is.null(bh.args$aa.mass)) {
      masses <- bh.args$aa.mass
      bh.args <- bh.args[ !('aa.mass' %in% names(bh.args)) ]
      if(!mass) {
        warning("incompatible arguments: forcing mass weighting")
        mass <- TRUE
      }
    }

    ## No mass-weighting
    else {
      masses <- NULL
      bh.args <- bh.args[ !('aa.mass' %in% names(bh.args)) ]
    }

    ## Build the Hessian Matrix
    if(is.null(hessian)) {
      cat(" Building Hessian...")
      ptm <- proc.time()
      H <- do.call('build.hessian', c(list(xyz=xyz, pfc.fun=ff, sequ=sequ, aa.mass=masses), bh.args))
      t <- proc.time() - ptm
      cat("\t\tDone in", t[[3]], "seconds.\n")
    }

    else {
      H <- hessian
    }


    ## Diagonalize matrix
    cat(" Diagonalizing Hessian...")
    ptm <- proc.time()
    ei <- eigen(H, symmetric=TRUE)
    t <- proc.time() - ptm
    cat("\tDone in", t[[3]], "seconds.\n")

    if(!is.null(keep)) {
      if(keep>ncol(ei$vectors))
        keep <- ncol(ei$vectors)
      keep <- keep-1
      keep.inds <- seq(ncol(ei$vectors)-keep, ncol(ei$vectors))
      ei$vectors <- ei$vectors[,keep.inds]
      ei$values <- ei$values[keep.inds]
    }

    ## Raw eigenvalues
    ei$values <- round(ei$values,6)
    triv.modes <- which(ei$values<=0) ## indicies !!

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
              This is useually an indication of an unphysical input structure.")
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
                call=cl)

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
