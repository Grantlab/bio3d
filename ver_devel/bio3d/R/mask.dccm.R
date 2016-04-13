mask <- function(...)
  UseMethod("mask")

mask.dccm <- function(dccm, pdb=NULL, a.inds=NULL, b.inds=NULL, ...) {
    x <- dccm
    dims <- dim(x)

    if(is.null(a.inds))
        stop("'a.inds' must be provided")

    if(is.null(b.inds))
        b.inds <- a.inds

    if(!is.null(pdb)) {
        if(!is.pdb(pdb))
            stop("If provided, 'pdb' must be a pdb object as obtained from 'read.pdb'")

        if(nrow(pdb$atom) != dims[1])
            stop("If provided, 'pdb' must match the dimensions of 'x'")

        if(!all(pdb$atom$elety == "CA"))
            warning("Non all-CA pdb detected. Ensure that provided 'pdb' match 'x'")

        if(is.select(a.inds))
            a.inds = a.inds$atom

        if(is.select(b.inds))
            b.inds = b.inds$atom

    }
    else {
        if(is.select(a.inds) | is.select(b.inds))
            stop("a.inds/b.inds should only be a 'select' object(s) when 'pdb' is provided")
    }

    if(!is.null(a.inds)) {
        tmp <- x
        tmp[ a.inds, b.inds ] = 666
        tmp[ b.inds, a.inds ] = 666
        inds.rm <- which(tmp != 666)

        x[inds.rm] = 0
    }

    return(x)
}
