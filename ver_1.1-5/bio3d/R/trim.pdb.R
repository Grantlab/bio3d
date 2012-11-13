`trim.pdb` <-
function(pdb, inds=NULL) {
  ## pdb <- read.pdb("http://www.rcsb.org/pdb/files/1poo.pdb")
  ## i <- atom.select(pdb,"back")
  ## n <- trim.pdb(pdb, inds=i)
  ## write.pdb(n, file="back.pdb")

  if(is.null(inds))
    stop("Selection indices, 'inds', from 'atom.select' required")

  if(!is.list(pdb))
    stop("Input 'pdb' must be a list object as returned from 'read.pdb'")

  new.pdb <- NULL
  new.pdb$atom <- pdb$atom[inds$atom,]
  new.pdb$xyz <-  pdb$xyz[inds$xyz]
  ##new.pdb$calpha <- as.logical(new.pdb$atom[, "elety"] == "CA")

  return(new.pdb)
}

