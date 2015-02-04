cmap.pdb <- function(pdb, inds=NULL, verbose=FALSE, ...) {
  if(!is.pdb(pdb))
    stop("provide a pdb object as obtained from function 'pdb'")

  if(is.null(inds)) {
    inds <- atom.select(pdb, "water", inverse=TRUE, verbose=verbose)
  }

  pdb <- trim.pdb(pdb, inds)
  xyz <- pdb$xyz
  grpby <- paste(pdb$atom$resno, pdb$atom$chain, sep="-")
  return(cmap.xyz(xyz, grpby, ...))
  

}
