#### Use Sys.setenv("PKG_CXXFLAGS"="-std=c++11")
#### when compiling

## maxlines=-1 missing
## ATOM.only=FALSE missing

read.pdb <- function(file, multi=FALSE, rm.insert=FALSE, rm.alt=TRUE,
                     ATOM.only = FALSE, verbose=FALSE, hex=FALSE) {
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
  pdb <- .read_pdb(file, multi=multi, hex=hex)
  if(!is.null(pdb$error))
    stop(paste("Could not read", file))
  else
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

  ## Remove 'Alt records'
  if (rm.alt) {
    if ( sum( !is.na(pdb$atom$alt) ) > 0 ) {
      first.alt <- sort( unique(na.omit(pdb$atom$alt)) )[1]
      cat(paste("   PDB has ALT records, taking",first.alt,"only, rm.alt=TRUE\n"))
      alt.inds <- which( (pdb$atom$alt != first.alt) ) # take first alt only
      if(length(alt.inds)>0) {
        pdb$atom <- pdb$atom[-alt.inds,]
        pdb$xyz <- trim.xyz(pdb$xyz, col.inds=-atom2xyz(alt.inds))
      }
    }
  }

  ## Remove 'Insert records'
  if (rm.insert) {
    if ( sum( !is.na(pdb$atom$insert) ) > 0 ) {
      cat("   PDB has INSERT records, removing, rm.insert=TRUE\n")
      insert.inds <- which(!is.na(pdb$atom$insert)) # rm insert positions
      pdb$atom <- pdb$atom[-insert.inds,]
      pdb$xyz <- trim.xyz(pdb$xyz, col.inds=-atom2xyz(insert.inds))
    }
  }
  
  if(any(duplicated(pdb$atom$eleno)))
    warning("duplicated element numbers ('eleno') detected")

  ## finished
  return(pdb)
}
