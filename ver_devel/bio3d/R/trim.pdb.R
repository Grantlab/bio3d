"trim.pdb" <-
  function(pdb, inds=NULL) {
    if(class(pdb)!="pdb")
      stop("input 'pdb' must be a list object as returned from 'read.pdb'")

    if(is.null(inds))
      stop("no selection indices provided")

    if(!is.list(inds))
      stop("selection indices must be provided i.e. from 'atom.select'")

    if(is.null(inds$atom) || is.null(inds$xyz))
      stop("selection indices must be provided i.e. from 'atom.select'")
        
    atom <- pdb$atom[inds$atom,]
    xyz <-  pdb$xyz[inds$xyz]
    calpha <- as.logical(atom[,"elety"]=="CA")
    
    if(!is.null(pdb$xyz.models))
      xyz.models <- pdb$xyz.models[,inds$xyz]
    else
      xyz.models <- NULL
    
    output<-list(atom=atom,
                 het=pdb$het, ## return unmodified
                 helix=pdb$helix, ## return unmodified
                 sheet=pdb$sheet, ## return unmodified
                 seqres=pdb$seqres, ## return unmodified
                 xyz=xyz,
                 xyz.models=xyz.models,
                 calpha = calpha)
    
    class(output) <- "pdb"
    return(output)
  }
