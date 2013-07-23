
"com" <-
  function(pdb, inds=NULL, use.mass=TRUE,
           mass.custom=NULL, elety.custom=NULL, rescue=TRUE) {
    
    if (missing(pdb)) 
      stop("com: must supply 'pdb' object, i.e. from 'read.pdb'")
    if(class(pdb)!="pdb")
      stop("com: 'pdb' must be of type 'pdb'")
    
    "atom2ele" <-
      function(atom.name, elety.dict, rescue=TRUE) {
        
        if( substr(atom.name,1,1) == "H" )
          atom.name <- "H"
        
        ele <- elety.dict[[ atom.name ]]
        
        if(is.null(ele)) {
          if(rescue) {
            ele <- substr(atom.name,1,1)
            warning(paste("unknown element: mapped ", atom.name, " to ", ele, sep=""))
          }
          else {
            stop(paste("atom2ele: element of '", atom.name, "' unknown", sep=""))
          }
        }
        return( ele )
      }
    
    "atom2mass" <-
      function(atom.name, mass.dict, elety.dict, rescue=TRUE) {
        
        ele <- atom2ele(atom.name, elety.dict, rescue=rescue)
        ele.mass <- mass.dict[[ ele ]]
        
        if( is.null(ele.mass) )
          stop(paste("atom2mass: mass of element '", ele, "' unknown", sep=""))
        return( ele.mass )
      }
    
    if(is.null(inds)) {
      xyz <- pdb$xyz
      at <- pdb$atom[, "elety"]
    }
    else {
      xyz <- pdb$xyz[inds$xyz]
      at <- pdb$atom[inds$atom, "elety"]
    }
    
    if(use.mass) {
      # Use LazyData to import data - changed Jul 23, 2013
#      data(atom.index)
      if(!is.null(elety.custom)) {
        if(class(elety.custom)!="list")
          stop("com: elety.custom must be of class 'list'")
        atom.index$elety <- c(elety.custom, atom.index$elety)
      }
      
      if (!is.null(mass.custom)) {
        if(class(mass.custom)!="list")
          stop("com: mass.custom must be of class 'list'")
        atom.index$mass <- c(mass.custom, atom.index$mass)
      }
      
      m <- unlist(lapply(at, atom2mass,
                         mass.dict = atom.index$mass,
                         elety.dict = atom.index$elety,
                         rescue=rescue))
    }
    else {
      m <- NULL
    }
    
    com <- com.xyz(xyz, m)
    out <- list("mass" = m, "xyz" = com)
    return(out)
  }
