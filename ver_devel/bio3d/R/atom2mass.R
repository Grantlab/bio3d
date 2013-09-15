"atom2mass" <-
  function(pdb, inds=NULL, mass.custom=NULL, elety.custom=NULL,
           grpby=NULL, rescue=TRUE) {
    
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
    
    if (!is.null(mass.custom)) {
      if(class(mass.custom)!="list")
        stop("mass.custom must be of class 'list'")
      atom.index$mass <- c(mass.custom, atom.index$mass)
    }
    
    ## fetch a list of element types
    eles <- atom2ele(atom.names, elety.custom=elety.custom, rescue=rescue)
    
    ele.mass <- NULL
    for ( ele in eles ) {
      ele.mass <- c(ele.mass, atom.index$mass[[ ele ]])
      
      if( is.null(ele.mass) )
        stop(paste("atom2mass: mass of element '", ele, "' unknown", sep=""))
    }
    
    if(!is.null(grpby)) {
      if( length(atom.names) != length(grpby) )
        stop("dimension miss-match in 'pdb' and 'grpby', check lengths")
      
      inds <- bounds(grpby, dup.inds=TRUE)
      ele.mass <- apply(inds, 1, function(ind, masses)
                        sum(masses[ind["start"]:ind["end"]]), ele.mass)
    }
    
    return( ele.mass )
  }
