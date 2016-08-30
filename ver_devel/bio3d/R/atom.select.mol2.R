.match.statbit <- function(pdb, statbit) {
  if(!is.character(statbit))
    stop("'statbit' must be a character vector")
  pdb$atom$statbit %in% statbit
}



atom.select.mol2 <- function(mol, string=NULL,
                             eleno = NULL, elety = NULL,
                             resid = NULL, chain = NULL, resno = NULL,
                             statbit = NULL, 
			     operator = "AND", inverse = FALSE,
                             value = FALSE, verbose=FALSE,  ...) {              

        
    #str.allowed <- c("noh", "h")
    #if(!is.null(string)) {
    #    if(!(string %in% str.allowed)) {
    #        stop("Not a valid selection string shortcut.\n\t Please use 'h' or 'noh'\n")
    #    }
    #}
    
    cl <- match.call()

    pdb <- as.pdb(mol)

    if(!is.mol2(mol))
        stop("'mol' must be an object of class 'mol2'")
    
    ## check input operator
    op.tbl <- c(rep("AND",3), rep("OR",4))
    operator <- op.tbl[match(operator, c("AND","and","&","OR","or","|","+"))]
    if(!operator %in% c("AND", "OR"))
        stop("Allowed values for 'operator' are 'AND' or 'OR'")
    
    ## check input string
    if(!is.null(string)) {
        str.allowed <- c("all", "protein", "notprotein", "nucleic", "notnucleic", "water", "notwater",
                         "calpha", "cbeta", "backbone", "back", "ligand", "h", "noh")
        if(!(string %in% str.allowed))
            stop("Unknown 'string' keyword. See documentation for allowed values")
    }

    ## verbose message output
    if(verbose) cat("\n")
    .verboseout <- function(M, type) {
        cat(" .. ", sprintf("%08s", length(which(M))), " atom(s) from '", type, "' selection \n", sep="")
    }
    
    ## combine logical vectors
    .combinelv <- function(L, M, operator) {
        if(operator=="AND") M <- L & M
        if(operator=="OR") M <- L | M
        return(M)
    }
    
    if(operator=="AND")
        M <- rep(TRUE, nrow(pdb$atom))
    if(operator=="OR")
        M <- rep(FALSE, nrow(pdb$atom))
    
    if(!is.null(string)) {   
        M <- switch(string,
                    all         =   M <- rep(TRUE, nrow(pdb$atom)),
                    protein     =  .is.protein(pdb),
                    notprotein  = !.is.protein(pdb),
                    nucleic     =  .is.nucleic(pdb),
                    notnucleic  = !.is.nucleic(pdb),
                    water       =  .is.water(pdb),
                    notwater    = !.is.water(pdb),
                    calpha      =  .is.protein(pdb)  & .match.elety(pdb, "CA"),
                    cbeta       =  .is.protein(pdb)  & .match.elety(pdb, c("CA", "N", "C", "O", "CB")),
                    backbone    =  .is.protein(pdb)  & .match.elety(pdb, c("CA", "N", "C", "O")),
                    back        =  .is.protein(pdb)  & .match.elety(pdb, c("CA", "N", "C", "O")),
                    ligand      = !.is.protein(pdb)  & !.is.nucleic(pdb) & !.is.water(pdb),
                    h           =  .is.hydrogen(pdb),
                    noh         = !.is.hydrogen(pdb),
                    NA
                    )
        
        if(verbose) {
            .verboseout(M, 'string')
        }
    }

  if(!is.null(eleno)) {
    L <- .match.eleno(mol, eleno)
    if(verbose) .verboseout(L, 'eleno')
    M <- .combinelv(L, M, operator)
  }
  if(!is.null(elety)) {
    L <- .match.elety(mol, elety)
    if(verbose) .verboseout(L, 'elety')
    M <- .combinelv(L, M, operator)
  }
  if(!is.null(resid)) {
    L <- .match.resid(mol, resid)
    if(verbose) .verboseout(L, 'resid')
    M <- .combinelv(L, M, operator)
  }
  if(!is.null(chain)) {
    L <- .match.chain(pdb, chain)
    if(verbose) .verboseout(L, 'chain')
    M <- .combinelv(L, M, operator)
  }
  if(!is.null(resno)) {
    L <- .match.resno(mol, resno)
    if(verbose) .verboseout(L, 'resno')
    M <- .combinelv(L, M, operator)
  }
  if(!is.null(statbit)) {
    L <- .match.statbit(mol, statbit)
    if(verbose) .verboseout(L, 'statbit')
    M <- .combinelv(L, M, operator)
  }

  if(verbose)
    cat(" ..", sprintf("%08s", length(which(M))), "atom(s) in final combined selection \n")

  if(inverse) {
    if(verbose) {
      cat(" ..", sprintf("%08s", length(which(!M))), "atom(s) in inversed selection \n")
    }
    sele <- as.select(which(!M))
  }
  else
    sele <- as.select(which(M))
    
  sele$call <- cl

  keep.bonds <- matrix(as.numeric(t(mol$bond[, c("origin","target")])) %in% sele$atom, ncol=2, byrow=T)
  bond.inds  <- which(apply(keep.bonds, 1, all))
  sele$bond <- bond.inds

  if(value)
    return(trim.pdb(pdb, sele))
  else
    return(sele)
}
