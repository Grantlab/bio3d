"atom.select.prmtop" <-
  function(prmtop, ...) {
    
    natoms <- prmtop$POINTERS[1]
    crds <- rep(NA, natoms*3)
    tmp.pdb <- amb2pdb(prmtop, crds)
    return(atom.select(tmp.pdb, ...))
}
