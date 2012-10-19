`pdbaln` <-
function(files, pqr=FALSE, ...) {

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

  missing <- !file.exists(files)
  if(any(missing)) {
    stop(paste(" ** Missing files:\n",
               paste( files[c(missing)], collapse="\n"),"\n",sep="") )
  }
  cat("Reading PDB files:")
  pdb.list <- NULL
  for(i in 1:length(files)) {
    if(pqr) {
      pdb.list[[ i ]] <- read.pqr(files[i])      
    } else {
      pdb.list[[ i ]] <- read.pdb(files[i])
    }
    cat(".")
  }
  
  cat("\n\nExtracting sequences\n")
  
  s <- lapply(pdb.list, seq.pdb)
  s <- t(sapply(s, `[`, 1:max(sapply(s, length))))
  s[is.na(s)] <- "-"
  s <- seqaln(s, id=files, ...)
  s <- read.fasta.pdb(s, pdb.path = "", pdbext = "")
  return(s)
}

