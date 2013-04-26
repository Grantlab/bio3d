`pdbaln` <-
function(files, pqr=FALSE, ncore=1, nseg.scale=1, ...) {

  ##- Quick and dirty alignment of pdb sequences
  ##   pdbs <- pdbaln(files)
  ##
  ##
  ## 'files' should be a character vector of input PDB file names
  ## '...' extra arguments for seqaln
  ##
  ## Improvements to include 'atom.select' arguments (chain
  ## spliting etc), formalisation of 'pdb.list' into a specific
  ## bio3d object of multiple structures like '3dalign'.
  ##
  ## pdb.list[[1]]$atom[1:3,]

  # Parallelized by multicore package (Fri Apr 26 19:24:18 EDT 2013)
  if(ncore > 1) {
     oops <- require(multicore)
     if(!oops)
        stop("Please install the multicore package from CRAN")

     options(cores = ncore)

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

  missing <- !file.exists(files)
  if(any(missing)) {
    stop(paste(" ** Missing files:\n",
               paste( files[c(missing)], collapse="\n"),"\n",sep="") )
  }
  cat("Reading PDB files:")
  mylapply <- lapply
  if(ncore > 1) mylapply <- mclapply
  pdb.list <- mylapply(1:length(files), function(i) {
    if(pqr) {
      pdb <- read.pqr(files[i])
    } else {
      pdb <- read.pdb(files[i])
    }
    cat(".")
    return( pdb )
  } )
  if(ncore > 1) readChildren()
#  pdb.list <- NULL
#  for(i in 1:length(files)) {
#    if(pqr) {
#      pdb.list[[ i ]] <- read.pqr(files[i])      
#    } else {
#      pdb.list[[ i ]] <- read.pdb(files[i])
#    }
#    cat(".")
#  }
  
  cat("\n\nExtracting sequences\n")
  
#  s <- lapply(pdb.list, seq.pdb)
  s <- mylapply(pdb.list, seq.pdb)
  if(ncore > 1) readChildren()

  s <- t(sapply(s, `[`, 1:max(sapply(s, length))))
  s[is.na(s)] <- "-"
  s <- seqaln(s, id=files, ...)
  s <- read.fasta.pdb(s, pdb.path = "", pdbext = "", ncore=ncore, nseg.scale=nseg.scale)
  return(s)
}

