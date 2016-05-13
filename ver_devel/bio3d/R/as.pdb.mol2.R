
as.pdb.mol2 <- function(mol2, ...) {
  
  natoms <- nrow(mol2$atom)
  xyz <- mol2$xyz
  tmp.pdb <- list()

  rownames(mol$substructure) <- mol$substructure$name
    resid <- mol$substructure[mol$atom$resid, "sub_type"]
  chain <- mol$substructure[mol$atom$resid, "chain"]

  tmp.pdb$atom <- data.frame(cbind(rep("ATOM", natoms),
                                   seq(1, natoms),
                                   mol2$atom$elena,
                                   NA,
                                   resid,
                                   chain, 
                                   mol2$atom$resno,
                                   NA,
                                   mol2$atom$x, mol2$atom$y, mol2$atom$z,
                                   NA, NA, NA,
                                   unlist(lapply(strsplit(mol2$atom$elety, split="[.]"),
                                                 function(x) x[1])),
                                   mol2$atom$charge),
                             stringsAsFactors=FALSE)

  
  colnames(tmp.pdb$atom) <- c("type", "eleno", "elety", "alt", "resid",
                              "chain", "resno", "insert",
                              "x", "y", "z", "o", "b", "segid", "elesy", "charge")
  
  tmp.pdb$xyz <- xyz
  class(tmp.pdb) <- "pdb"
  ##ca.inds        <- rep(FALSE, natoms)
  return(tmp.pdb)
}
