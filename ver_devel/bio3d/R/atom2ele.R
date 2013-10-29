"atom2ele" <-
  function(pdb, inds=NULL, elety.custom=NULL, rescue=TRUE) {

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
    
    if(!is.null(elety.custom)) {
      if(class(elety.custom)!="list")
        stop("elety.custom must be of class 'list'")
      atom.index$elety <- c(elety.custom, atom.index$elety)
    }

    atom.names[substr(atom.names,1,1) == "H"] <- "H"
    eles <- atom.index$elety[atom.names]
    is.unknown <- is.na(names(eles))

    if(any(is.unknown)) {
      if(rescue) {
        eles[is.unknown] <- substr(atom.names[is.unknown],1,1)
        warning(paste("\n\tunknown element: mapped ", atom.names[is.unknown], " to ", eles[is.unknown], sep=""))
      }
      else {
        stop(paste("\n\tatom2ele: element of '", atom.names[is.unknown], "' unknown", sep=""))
      }
    }
    eles <- unlist(eles)
    names(eles) <- NULL
    return( eles )
  }
