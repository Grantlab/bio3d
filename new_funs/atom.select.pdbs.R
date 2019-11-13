atom.select.pdbs <- function(pdbs,
                             string = NULL, 
                             resno = NULL, chain = NULL, resid = NULL,
                             operator="AND", inverse = FALSE,
                             value = FALSE, verbose=FALSE,  ...) {

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
    
    cl <- match.call()
    if(operator=="AND")
        M <- rep(TRUE, ncol(pdbs$ali))
    if(operator=="OR")
        M <- rep(FALSE, nrow(pdbs$ali))


    if(!is.null(string)) {
        gaps <- gap.inspect(pdbs$ali)
        M <- switch(string,
                    all         =  M <- rep(TRUE, ncol(pdbs$ali)),
                    gap         =  apply(pdbs$ali, 2, function(x) any(x %in% "-")),
                    gaps        =  apply(pdbs$ali, 2, function(x) any(x %in% "-")),
                    nongap      = !apply(pdbs$ali, 2, function(x) any(x %in% "-")),
                    allgap      =  apply(pdbs$ali, 2, function(x) all(x %in% "-"))
                    )

        if(verbose) {
            .verboseout(M, 'string')
        }
    }
    
    if(!is.null(resno)) {
        L <- apply(pdbs$resno, 2, function(x) any(x %in% resno))
        if(verbose) .verboseout(L, 'resno')
        M <- .combinelv(L, M, operator)
    }

    if(!is.null(chain)) {
        L <- apply(pdbs$chain, 2, function(x) any(x %in% chain))
        if(verbose) .verboseout(L, 'chain')
        M <- .combinelv(L, M, operator)
    }

    if(!is.null(resid)) {
        L <- apply(pdbs$resid, 2, function(x) any(x %in% resid))
        if(verbose) .verboseout(L, 'resid')
        M <- .combinelv(L, M, operator)
    }
    
    if(inverse) {
        if(verbose) {
            cat(" ..", sprintf("%08s", length(which(!M))), "atom(s) in inversed selection \n")
        }
        sele <- as.select(which(!M))
    }
    else
        sele <- as.select(which(M))
    
    sele$call <- cl
    if(verbose) cat("\n")
    
    if(value)
        return(trim(pdbs, col.inds=sele$atom))
    else
        return(sele)
}
