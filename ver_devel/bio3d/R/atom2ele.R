"atom2ele" <-
  function(pdb, inds=NULL, elety.custom=NULL, rescue=TRUE) {
    
    if(class(pdb)=="pdb") {
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
    
    eles <- NULL
    for ( at in atom.names ) {
      if( substr(at,1,1) == "H" )
        at <- "H"
      ele <- atom.index$elety[[ at ]]
      
      if(is.null(ele)) {
        if(rescue) {
          ele <- substr(at, 1,1)
          warning(paste("unknown element: mapped ", at, " to ", ele, sep=""))
        }
        else {
          stop(paste("atom2ele: element of '", at, "' unknown", sep=""))
        }
      }
      eles <- c(eles, ele)
    }
    return( eles )
  }
