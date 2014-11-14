"atom.select.prmtop" <-
  function(prmtop, ...) {
    
    tmp.pdb <- amb2pdb(prmtop, crd=NULL)
    return(atom.select.pdb(tmp.pdb, ...))
}
