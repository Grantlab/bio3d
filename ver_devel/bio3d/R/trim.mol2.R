trim.mol2 <- function(mol, inds=NULL, ...) {

    if(!is.mol2(mol))
        stop("input should be of class 'mol2' as obtained by 'read.mol2'")
    
    if(!is.select(inds))
        stop("provide indices")

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
    new$substructure <- mol$substructure[ mol$substructure$root_atom %in% mol$atom$eleno[ inds$atom], ]
    new$substructure$root_atom <- new$atom[ as.character(new$substructure$root_atom), "eleno" ]

    ## new ids and rownames
    new$bond$id <- 1:nrow(new$bond)
    new$substructure$id <- 1:nrow(new$substructure)
    
    rownames(new$bond) <- 1:nrow(new$bond)
    rownames(new$atom) <- 1:nrow(new$atom)
    rownames(new$substructure) <- 1:nrow(new$substructure)
  
    return(new)
 }

    
