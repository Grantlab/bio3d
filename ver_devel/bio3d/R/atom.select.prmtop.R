"atom.select.prmtop" <-
  function(prmtop, ...) {
    
    
    tmp.pdb <- amb2pdb(prmtop, crds=NULL)
    return(atom.select.pdb(tmp.pdb, ...))
}
