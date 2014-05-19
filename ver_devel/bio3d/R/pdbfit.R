pdbfit <-
function(pdbs, inds=NULL, outpath=NULL, ...) {
  ##
  ## Quick Fit Fitter for PDBs
  ##  was called 'fit.pdbs()' in model.R
  ##
  if(!inherits(pdbs, "3dalign")) {
    stop("Input 'pdbs' should be of class '3dalign', e.g. from pdbaln() or read.fasta.pdb()")
  }
  full <- ifelse(is.null(outpath), FALSE, TRUE)
  if(is.null(inds)) {    inds <-gap.inspect(pdbs$xyz)$f.inds  }
  if(is.list(inds)){ inds=inds$xyz }
  return( fit.xyz( fixed=pdbs$xyz[1,], mobile=pdbs, fixed.inds=inds,
                  mobile.inds=inds, outpath=outpath,
                  full.pdbs=full, ... ))
}
