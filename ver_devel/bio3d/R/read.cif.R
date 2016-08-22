read.cif <- function(file, maxlines = -1, multi=FALSE, rm.insert=FALSE, rm.alt=TRUE,
                     ATOM.only=FALSE, verbose=FALSE) {

  cl <- match.call()
  warning("helix/sheet records will not be parsed in this version of the code")

  if(missing(file)) {
    stop("read.cif: please specify a CIF 'file' for reading")
  }
  
  if(!is.logical(multi)) {
    stop("read.cif: 'multi' must be logical TRUE/FALSE")
  }

  ##- Check if file exists locally or on-line
  toread <- file.exists(file)
  if(substr(file,1,4)=="http") {
    toread <- TRUE
  }
      
  ## Check for 4 letter code and possible on-line file
  if(!toread) {
    if(nchar(file)==4) {
      cat("  Note: Accessing on-line CIF file\n")
      file <- get.pdb(file, path=tempdir(), quiet=TRUE, format="cif")
    }
    else {
      stop("No input CIF file found: check filename")
    }
  }
  else {
      file <- normalizePath(file)
  }
  
  ## parse CIF file with cpp function
  pdb <- .read_cif(file, multi=FALSE, hex=FALSE, maxlines=maxlines, atoms_only=ATOM.only)
  if(!is.null(pdb$error))
    stop(paste("Error in reading CIF file", file))
  else
    class(pdb) <- c("pdb") ## this should be cif perhaps?

  if(pdb$models > 1)
      pdb$xyz <- matrix(pdb$xyz, nrow=pdb$models, byrow=TRUE)
  
  pdb$xyz <- as.xyz(pdb$xyz)

  ## set empty strings to NA
  pdb$atom[pdb$atom==""] <- NA
  pdb$atom[pdb$atom=="?"] <- NA
  pdb$atom[pdb$atom=="."] <- NA

  ## set call
  pdb$call <- cl
  
  return(pdb)
}

