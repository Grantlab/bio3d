"atom2mass" <-
  function(pdb, inds=NULL, mass.custom=NULL, elety.custom=NULL,
           grpby=NULL, rescue=TRUE) {
    
    if(is.pdb(pdb)) {
      if(!is.null(inds)) {
        pdb <- trim.pdb(pdb, inds)
      }
      atom.names <- pdb$atom[,"elety"]
    }
    else {
      if(!is.null(inds))
        warning("'inds' not in use when 'pdb' is vector")
      atom.names <- pdb
    }
    
    if (!is.null(mass.custom)) {
      if(class(mass.custom)!="list")
        stop("mass.custom must be of class 'list'")
      atom.index$mass <- c(mass.custom, atom.index$mass)
    }
    
    ## fetch a list of element types
    eles <- atom2ele(atom.names, elety.custom=elety.custom, rescue=rescue)

    ele.mass <- unlist(atom.index$mass[eles])
    
    if(any(is.na(ele.mass)))
        stop(paste("\n\tatom2mass: mass of element '", eles[is.na(ele.mass)], "' unknown", sep=""))
    
    if(!is.null(grpby)) {
      if( length(atom.names) != length(grpby) )
        stop("dimension miss-match in 'pdb' and 'grpby', check lengths")
      ele.mass <- unlist(lapply(split(ele.mass, grpby), sum))
    }
    names(ele.mass) <- NULL
    return( ele.mass )
  }
