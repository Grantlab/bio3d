pdbaln <-
function(files, fit=FALSE, pqr=FALSE, ncore=1, nseg.scale=1, progress=NULL, ...) {

  ## Log the call
  cl <- match.call()

  ##- Quick and dirty alignment of pdb sequences
  ##   pdbs <- pdbaln(files)
  ##
  ##
  ## 'files' should be a character vector of input PDB file names
  ## '...' extra arguments for seqaln
  ##
  ## Improvements to include 'atom.select' arguments (chain
  ## spliting etc), formalisation of 'pdb.list' into a specific
  ## bio3d object of multiple structures like 'pdbs'.
  ##
  ## pdb.list[[1]]$atom[1:3,]

  # Parallelized by parallel package (Fri Apr 26 19:24:18 EDT 2013)
  ncore <- setup.ncore(ncore)
  if(ncore > 1) {
     # Issue of serialization problem
     # Maximal number of cells of a double-precision matrix
     # that each core can serialize: (2^31-1-61)/8
     R_NCELL_LIMIT_CORE = 2.68435448e8
     R_NCELL_LIMIT = ncore * R_NCELL_LIMIT_CORE

     if(nseg.scale < 1) {
        warning("nseg.scale should be 1 or a larger integer\n")
        nseg.scale=1
     }
  }
  
  if(!is.list(files)) {
    ## Check if input PDB files exist localy or online
    toread.local <- file.exists(files)
    toread.online <- (substr(files,1,4)=="http")
    toread.id <- rep(FALSE, length(files))
    toread <- as.logical(toread.local + toread.online)
    
    ## Check for 4 letter code and possible online file
    if(any(!toread)) {
      toread.id <- ((nchar(files)==4) + (!toread) == 2)
      files[toread.id] <-  get.pdb(files[toread.id], URLonly=TRUE)
      cat("  Note: Accessing online PDB files using 4 letter PDBID\n")
    }
    
    ## Exit if we still have missing files
    missing <- !as.logical(toread + toread.id)
    if(any(missing)) {
      stop(paste(" ** Missing files: check filenames\n",
                 paste( files[c(missing)], collapse="\n"),"\n",sep="") )
    }
    
    ## Avoid multi-thread downloading
    if(any(toread.online | toread.id)) {
      ncore = 1
#      options(cores = ncore)
    }
    cat("Reading PDB files:",files, sep="\n")
    
    pdb.list <- mclapply(1:length(files), function(i) {
      if(pqr) {
        pdb <- read.pqr(files[i])
      } else {
        pdb <- read.pdb(files[i])
      }
      cat(".")
      
      if(!is.null(progress)) {
        progress$inc(1/length(files)/2)
      }

      return( pdb )
    } , mc.cores = ncore)
  }
  else {
    if(!all(sapply(files, is.pdb)))
      stop("'files' must be a vector of file names, or a list of pdb objects")
    pdb.list <- files
    #ids <- sapply(files, function(x) x$call$file)
    #print(ids)
  }

  cat("\n\nExtracting sequences\n")
  s <- mclapply(pdb.list, pdbseq, mc.cores=ncore)

  ## check for NULL in pdbseq output
  ## (this would indicate no amino acid sequence in PDB)
  tmpcheck <- unlist(lapply(s, is.null))
  if(any(tmpcheck)) {
    if(!is.list(files)) {
      err <- paste(
        "Could not align PDBs due to missing amino acid sequence in files:\n ",
        paste(files[tmpcheck], collapse=", ")
        )
    }
    else {
      err <- paste(
        "Could not align PDBs due to missing amino acid sequence in PDB:\n ",
        paste(tmpcheck, collapse=", ")
        )
    }
    stop(err)
  }

  s <- t(sapply(s, `[`, 1:max(sapply(s, length))))
  s[is.na(s)] <- "-"
  
  if(!is.list(files)) {
    s <- seqaln(s, id=files, ...)
    files <- NULL
  }
  else {
    s <- seqaln(s, id=NULL, ...)
  }

  cat("\n")
  s <- read.fasta.pdb(s, pdblist=files,
                      prefix = "", pdbext = "",
                      ncore=ncore, nseg.scale=nseg.scale, progress=progress)
  s$call=cl
  
  if(fit)
    s$xyz <- pdbfit(s)
  return(s)
}

