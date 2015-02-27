#### Use Sys.setenv("PKG_CXXFLAGS"="-std=c++11")
#### when compiling

read.pdb2 <- function(file) {
  pdb <- .read_pdb(file)
  
  if(pdb$models>1)
    pdb$xyz <- matrix(pdb$xyz, nrow=pdb$models, byrow=TRUE)
  
  pdb$xyz <- as.xyz(pdb$xyz)
  pdb$models <- NULL
  class(pdb) <- "pdb"
  
  ca.inds <-  atom.select.pdb(pdb, string="calpha", verbose=FALSE)
  pdb$calpha <- seq(1, nrow(pdb$atom)) %in% ca.inds$atom
  pdb$atom[pdb$atom==""] <- NA

  warning("ignoring SSE and SEQRES records in PDB file")
  return(pdb)
}
