"binding.site" <-
  function(a, b = NULL,a.inds = NULL, b.inds = NULL,
           cut=5, hydrogens=TRUE) {

  if (missing(a))
    stop("binding.site: must supply an input 'pdb' object 'a', i.e. from 'read.pdb'")

  if(!is.null(b)) {

    if ( hydrogens ) {
      if(is.null(a.inds))
        a.inds <- atom.select(a, string='///////')
      if(is.null(b.inds))
        b.inds <- atom.select(b, string='///////')
    }
    else {
      if(is.null(a.inds))
        a.inds <- atom.select(a, string='noh')
      if(is.null(b.inds))
        b.inds <- atom.select(b, string='noh')
    }
  }

  else {
    complex <- a
    a <- trim.pdb(complex, a.inds)
    b <- trim.pdb(complex, b.inds)

    if ( hydrogens ) {
      a.inds <- atom.select(a, string='///////')
      b.inds <- atom.select(b, string='///////')
    }
    else {
       a.inds <- atom.select(a, string='noh')
       b.inds <- atom.select(b, string='noh')
     }
  }

  # Join the coordinates of the two entities
  c <- c(a$xyz[a.inds$xyz], b$xyz[b.inds$xyz])
  last <- as.numeric(a$atom[nrow(a$atom),"resno"])

  # .. and make the last PDB entity to one residue
  b$atom[,"resno"] <- last+1

  # Calcualte distance matrix and group by residue number
  dmat <- dm.xyz(c, grpby=c(a$atom[a.inds$atom,"resno"],
                      b$atom[b.inds$atom,"resno"]), scut=0)

  resno.map <- unique(c(a$atom[a.inds$atom,"resno"], b$atom[b.inds$atom,"resno"]))
  ligresno <- unique(b$atom[b.inds$atom,"resno"])

  # Fetch ligand residue indices and its distances to the protein
  inds <- which(resno.map %in% ligresno)
  distances <- dmat[,inds]
  close.inds <- which(distances<cut)

  # Make the output
  atom.inds <- which(a$atom[,"resno"] %in% resno.map[close.inds])
  atom.inds <- intersect(atom.inds, a.inds$atom)
  xyz.inds <- atom2xyz(atom.inds)

  resno <- as.numeric( unique(a$atom[atom.inds, "resno"]) )
  resnames <- apply(a$atom[atom.inds,c("resid", "resno")], 1, paste, collapse="")
  resnames <- unique(resnames)


  out <- list("atom.inds"=atom.inds, "xyz.inds"=xyz.inds,
              "resnames"=resnames, "resno"=resno)
  return(out)
}

