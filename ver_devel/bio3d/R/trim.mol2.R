trim.mol2 <- function(mol, ..., inds=NULL) {

    if(!is.mol2(mol))
        stop("input should be of class 'mol2' as obtained by 'read.mol2'")
    
    cl <- match.call()
    extra.args <- list(...)
    
    if(length(extra.args)>0) {
      if(!is.null(inds))
          warning("Multiple atom selection terms provided. Using only argument 'inds'")
      else if(is.select(extra.args[[1]]))
          ## to be back-compatible with the habit calling trim.pdb(pdb, inds)
          inds = extra.args[[1]]
      else
          inds = atom.select(mol, ...)
    }

    if(is.null(inds))
        stop("no selection indices provided")
    
    if(!is.list(inds))
      stop("selection indices must be provided i.e. from 'atom.select'")
    
    if(is.null(inds$atom) || is.null(inds$xyz))
        stop("selection indices must be provided i.e. from 'atom.select'")
    

    new <- mol
    rownames(new$atom) <- new$atom$eleno
    new$atom <- new$atom[inds$atom,,drop=FALSE]
    new$atom$eleno <- 1:nrow(new$atom)
    new$bond <- mol$bond[inds$bond, ]
    
    xyz <- as.xyz(mol$xyz)
    new$xyz <- as.xyz(mol$xyz[, inds$xyz, drop=FALSE])
    
    new$info[1] <- nrow(new$atom)
    new$info[2] <- nrow(new$bond)
    new$info[3] <- length(unique(new$atom$resid))

    new$bond$origin = new$atom[ as.character(new$bond$origin), "eleno" ]
    new$bond$target = new$atom[ as.character(new$bond$target), "eleno" ]

    ## substrucutre
    if(!is.null(mol$substructure)) {
        new$substructure <- mol$substructure[ mol$substructure$root_atom %in% mol$atom$eleno[ inds$atom], ]
        
        if(nrow(new$substructure) > 0) {
            new$substructure$root_atom <- new$atom[ as.character(new$substructure$root_atom), "eleno" ]
            new$substructure$id <- 1:nrow(new$substructure)
        }
    }

    ## new ids and rownames
    new$bond$id <- 1:nrow(new$bond)
    rownames(new$bond) <- 1:nrow(new$bond)
    rownames(new$atom) <- 1:nrow(new$atom)
 
    return(new)
 }

    
