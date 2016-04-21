

read.hessian <- function(filename) {
    hess <- as.numeric( readLines(filename) )
    natom <- sqrt(length(hess) / 9)
    h <- matrix(hess, ncol=3*natom, nrow=3*natom, byrow=T)
    return(h)
}

effective.hessian <- function(h, inds) {
    haa <- h[inds$xyz, inds$xyz]
    hao <- h[inds$xyz, -inds$xyz]
    hoo <- h[-inds$xyz, -inds$xyz]
    hoa <- t(hao)

    hoo.inv <- solve(hoo)
    k <- haa - ((hao %*% hoo.inv) %*% hoa)
    return(k)
}

build.modes <- function(k, xyz) {

    ei <- eigen(k)
    U <- ei$vectors[, seq(ncol(ei$vectors), 1)]
    L <- rev(ei$values)
    head(L, n=12)

    modes <- list(modes = U,
                  U = U,
                  L = L,
                  frequencies=NULL,
                  force.constants=L,
                  fluctuations=NULL,
                  temp=NULL, masses=NULL,
                  xyz = xyz,
                  triv.modes = 6,
                  natoms = length(xyz)/3,
                  call = NULL)

    class(modes) <- c("EnergeticModes", "nma")
    modes$fluctuations <- fluct.nma(modes, mode.inds=NULL)
    return(modes)
}

hessian2fc <- function(h, mask.upper=TRUE) {
    dims <- dim(h)
    natoms <- dims[1]/3
    k.mat <- matrix(0, ncol=natoms, nrow=natoms)

    rinds <- seq(1, natoms)
    for ( i in 1:natoms ) {
        inds.a <- atom2xyz(i)
        for ( j in 1:natoms ) {
            inds.b <- atom2xyz(j)
            tmp <- h[inds.a, inds.b]
            kij <- sum(diag(tmp))

            k.mat[i,j] <- kij
        }
    }

    k.mat <- (-1) * k.mat
    kmat.tri <- k.mat
    kmat.tri[ lower.tri(k.mat, diag=TRUE) ] <- NA

    return(kmat.tri)
}

build.dmat <- function(xyz) {
    r.mat <- dm.xyz(xyz, mask.lower=TRUE)
    diag(r.mat) <- NA
    return(r.mat)
}

mat2vec <- function(mat) {
    mat = c(t(mat))
    return(mat[ !is.na(mat) ])
}

atompairs <- function(pdb, inds=NULL) {
    if(is.null(inds))
        pdb <- trim(pdb, inds)

    labs <- paste(pdb$atom$resid, "-",
                  pdb$atom$resno, '@', pdb$atom$elety, sep='')

    plabs = c()
    N = nrow(pdb$atom)

    for(i in 1:(N-1)) {
        inds.b = (i+1):N
        inds.a = rep(i, length(inds.b))

        labs.a = labs[inds.a]
        labs.b = labs[inds.b]

        tmp = paste(labs.a, labs.b)
        plabs = c(plabs, tmp)
    }

    return(plabs)
}
