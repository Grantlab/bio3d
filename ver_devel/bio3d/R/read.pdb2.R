#### Use Sys.setenv("PKG_CXXFLAGS"="-std=c++11")
#### when compiling

read.pdb2 <- function(file, multi=FALSE) {
  cl <- match.call()

  if(missing(file)) {
    stop("read.pdb: please specify a PDB 'file' for reading")
  }
  
  if(!is.logical(multi)) {
    stop("read.pdb: 'multi' must be logical TRUE/FALSE")
  }

  ##- Check if file exists locally or on-line
  toread <- file.exists(file)
  if(substr(file,1,4)=="http") {
    toread <- TRUE
  }

  ## Check for 4 letter code and possible on-line file
  if(!toread) {
    if(nchar(file)==4) {
      cat("  Note: Accessing on-line PDB file\n")
      file <- get.pdb(file, path=tempdir(), quiet=TRUE)
    }
    else {
      stop("No input PDB file found: check filename")
    }
  }
  
  ## parse PDB file with cpp function
  pdb <- .read_pdb(file, multi=multi)
  class(pdb) <- "pdb"

  ## convert xyz to matrix
  if(pdb$models>1)
    pdb$xyz <- matrix(pdb$xyz, nrow=pdb$models, byrow=TRUE)
  
  pdb$xyz <- as.xyz(pdb$xyz)
  pdb$models <- NULL

  ## construct c-alpha attribute
  ca.inds <-  atom.select.pdb(pdb, string="calpha", verbose=FALSE)
  pdb$calpha <- seq(1, nrow(pdb$atom)) %in% ca.inds$atom

  ## set empty strings to NA
  pdb$atom[pdb$atom==""] <- NA

  ## give names to seqres
  names(pdb$seqres) <- pdb$seqres_chain
  pdb$seqres_chain <- NULL

  ## set call
  pdb$call <- cl

  ## finished
  return(pdb)
}
