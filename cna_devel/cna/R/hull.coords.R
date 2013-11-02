hull.coords <- function(membership,
                             pdb
                             ){

  ## Check if the pdb object has 'pdb' class
  oops <- class(pdb)
  if(sum(oops %in% "pdb")<0){
    stop("Error! pdb object must be a 'pdb' class object")
  }

  ## Do a renumbering of pdb resno to match comms$membership and indices
  pdb <- convert.pdb(pdb, type="pdb", renumber=TRUE)

  res.num <- length(pdb$atom[,"resno"])

  if(res.num != length(membership)){
    stop("The pdb file and the membership vector provided have different number of elements")
  }
  
  ## Select the atoms for each community and calculate its
  ## geometric center
  coords <- list()

  for(a in 1:max(membership)){
    residues <- which(membership==a)
    noh.inds <- atom.select(pdb, resno=residues, string="noh", verbose=FALSE)

    coords.table <- matrix(pdb$xyz[noh.inds$xyz],
                           nrow = dim(pdb$atom[noh.inds$atom,,drop=FALSE])[1],
                           ncol = 3,
                           byrow=TRUE)
 
    #coords[[a]] <- chull(coords.table)
    coords[[a]] <- convhull(coords.table)
  }

  return(coords)
  
}
           
