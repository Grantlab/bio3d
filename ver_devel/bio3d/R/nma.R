"nma" <-
  function(pdb, inds=NULL, ff='calpha', pfc.fun=NULL,
           normalize=TRUE, mass=TRUE, temp=300.0, keep=NULL,
           compiler=FALSE, ncore=1, cutoff=15, gamma=1  ) {
    
    if (missing(pdb))
      stop("nma: must supply 'pdb' object, i.e. from 'read.pdb'")
    if(class(pdb)!="pdb")
      stop("nma: 'pdb' must be of type 'pdb'")
    
    ## Log the call
    cl <- match.call()
    
    ## Trim PDB to match user selection
    if ( !is.null(inds) ) {
      pdb <- trim.pdb(pdb, inds)
      ## can be removed with new version of trim.pdb
      ##pdb$calpha <- as.logical(pdb$atom[,"elety"] == "CA")
    }
    
    ## Define force field
    if (is.null(pfc.fun)) {
      
      ## Bahar "ANM"-ff
      if (ff=="anm")  {
        if(normalize){
          warning("nma: set 'normalize=FALSE' when using force field 'anm'")
          ##normalize <- FALSE
        }
        "ff.anm" <- function(r, rc=cutoff, g=gamma) {
          ifelse( r>rc, 0, g/(r**2) )
        }
        ff <- ff.anm
      }
      
      ## Hinsen "C-alpha"-ff
      else if (ff=="calpha")  {
        if(!normalize){
          warning("nma: set 'normalize=TRUE' when using force field 'calpha'")
          ##normalize <- TRUE
        }
        "ff.calpha" <- function(r) {
          a <- 1e-1; b <- 1; c <- 1e6;
          ifelse( r<4.0,
                 (a*8.6*(10**5)*r) - (b*2.39*(10**5)), 
                 c*128 * r**(-6) )
        }
        ff <- ff.calpha
      }
      
      else
        stop("nma: options for 'ff' is 'calpha' or 'anm'")
      
    }
    else {
      ## Use customized force field
      ff <- pfc.fun
    }
    
    ## Define residues masses
    ## Source: MMTK (for reproduction purposes!)
    w <- c( 71.079018, 157.196106, 114.104059, 114.080689, 103.143407,
           128.107678, 128.131048,  57.05203,  137.141527, 113.159985,
           113.159985, 129.18266,  131.197384, 147.177144,  97.117044,
           87.078323, 101.105312, 186.213917, 163.176449,  99.132996)
    
    aa <- c("ALA", "ARG", "ASN", "ASP", "CYS",
            "GLU", "GLN", "GLY", "HIS", "ILE",
            "LEU", "LYS", "MET", "PHE", "PRO",
            "SER", "THR", "TRP", "TYR", "VAL")

    "res2wt" <- function(x, w, aa) {
      ind <- which(aa==x)
      if(length(ind)==1)
        return(w[ind])
      else
        return(NA)
    } 
    
    ## Only C-alpha ENM NMA is implemented
    ca.inds <- atom.select(pdb, "calpha", verbose=FALSE)
    natoms <- length(ca.inds$atom)
    if (natoms<3)
      stop("nma: insufficient number of CA atoms in structure")
    xyz <- pdb$xyz[ca.inds$xyz]
    
    ## If mass-weighted force constant matrix
    if (mass) {
      sequ <- pdb$atom[ca.inds$atom,"resid"]
      wts <- unlist(lapply(sequ, res2wt, w, aa))
      if(NA%in%wts)
        stop("nma: unknown residue type")
      wts <- sqrt(wts)
    } else {
      wts <- NULL
    }
    
    ## Build the Hessian Matrix - use byte compiled code if possible
    cat(" Building Hessian...")
    ptm <- proc.time()
    H <- build.hessian(xyz, ff, normalize=normalize, mass.weights=wts, 
                       compiler=compiler, ncore=ncore)
    t <- proc.time() - ptm
    cat("\t\tDone in", t[[3]], "seconds.\n")
  
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
    triv.modes <- which(ei$values==0) ## indicies !!
    
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
      tri.inds <- rep(1:natoms, each=3)
      V <- apply(V, 2, '*', 1 / wts[tri.inds])
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
                mass=wts,
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
