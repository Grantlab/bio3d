## perturbation matrix
## row i, avg displacement for residue i when all other residues  are perturbed one at the time
## col j, response profile of all residues under the perturbation of residue j
prs1 <- function(h.inv, fv=NULL, reps=10) {
    natoms <- nrow(h.inv) / 3
    pm <- matrix(0, nrow=natoms, ncol=natoms)

    ## calculate random force vector
    if(is.null(fv))
      fv <- forcegen.default(reps)

    for(i in 1:natoms) {

        for(j in 1:nrow(fv)) {
            
            ## apply 'kick'
            r <- h.inv[, atom2xyz(i)] %*% fv[j,]

            ## amplitudes for each residue
            d <- matrix(r, ncol=3, byrow=TRUE)**2
            d2 <- sqrt(rowSums(d))
            pm[,i] <- pm[,i] + d2
        }

    }

    pm <- pm / reps
    return(pm)

}

## same as prs1 but with grpby argument
## this function is ment to perturb all atoms within a residue
prs3 <- function(h.inv, fv=NULL, reps=10, grpby=NULL) {

  if(is.null(grpby)) {
    natoms <- nrow(h.inv) / 3
  }
  else {
    natoms <- length(unique(grpby))
    grpby.unq <- unique(grpby)
  }
  
  pm <- matrix(0, nrow=natoms, ncol=natoms)

  ## calculate random force vector
  if(is.null(fv))
    fv <- forcegen.default(reps)

  for(i in 1:natoms) {
    if(!is.null(grpby))
      k <- which(grpby == grpby.unq[i])
    else
      k <- i
    
    for(j in 1:nrow(fv)) {
      if(length(k) > 1)
        v <- rep(fv[j,], times=length(k))
      else
        v <- fv[j,]
      
      ## apply 'kick'
      r <- h.inv[, atom2xyz(k)] %*% v

      ## amplitudes for each residue
      d <- matrix(r, ncol=3, byrow=TRUE)**2
      d2 <- sqrt(rowSums(d))
      d3 <- unlist(lapply(grpby.unq, function(x) sum(d2[ grpby == x ])))
      
      pm[,i] <- pm[,i] + d3
    }

  }
  
  pm <- pm / reps
  return(pm)
  
}



## DFI(i): totalt displacement of residue i,
## induced by perturbations placed on all residues in the protein
dfi <- function(A) {
    d <- colSums(A) / sum(A)
    return(d)
}

## allosteric response ratio
arr <- function(A, inds) {
    d = rep(0, ncol(A))

    for(j in 1:ncol(A)) {

        d[j] = sum(A[,inds]) / length(inds)
        d[j] = d[j] / (sum(A[,j]) / (ncol(A)-1))
        
    }
    return(d)
    
}

## returns a matrix of N random vectors
forcegen.default <- function(reps=10) {

    ## generate N sets of random vectors
    fv <- matrix(NA, nrow=reps, ncol=3)
    for(i in 1:nrow(fv))
        fv[i,] <- runif(3, 0, 1)

    ## normalize force vectors
    F <- fv / sqrt(rowSums(fv**2))
    return(F)
    
}



## not use
forcegen.pdb <- function(pdb, sele) {
    ## random force vector (unit vector)
    v <- runif(3, 0, 1)
    v2 <- rep(v, times=length(sele$xyz)/3)
    fv <- v2 / sqrt(sum(v2**2))

    N <- nrow(pdb$atom)
    F <- rep(0, 3*N)
    F[ sele$xyz ] <- fv
    return(F)
}

## wrapper / utility function for plotting a standard fluctuation plot
## with colored regions ('regs') in one function call
plotb2 <- function(x, ..., regs=NULL,
                   col.reg=rep("lightblue", nrow(regs))) {

    ylim=range(x)
    xlim=c(1, length(x))


    plot.new()
    plot.window(xlim=xlim, ylim=ylim)

    rect(regs[,"start"], rep(ylim[1], nrow(regs)), regs[,2],
         rep(ylim[2], nrow(regs)),
         col=col.reg,
         border=NA)

    ## add this for plot.bio3d on the same device
    par(new=TRUE)

    plot.bio3d(x, ...)

}


## wrapper / utility function for plotting perturbation matrix
## with annotations (sse and regions) in one function call
image3 <- function(x=NULL, y=NULL, z=NULL,
                   col=colorRampPalette(c("white", "red"))(10),
                   xlim=NULL, ylim=NULL,
                   xlab="Resno", ylab=xlab,
                   resno = NULL,
                   sse = NULL, regs = NULL,
                   col.reg = "lightblue", ...) {

    if(is.null(x))
        x=1:nrow(z)
    if(is.null(y))
        y=1:ncol(z)

    if(is.null(resno))
        resno <- x
    
    if(is.null(xlim))
        xlim=range(x)
    if(is.null(ylim))
        ylim=range(y)

    xlim = xlim + c(-5, 5)
    ylim = ylim + c(-5, 5)


    image(x, y, z,
          xlim = xlim, ylim = ylim,
          xlab = xlab, ylab = ylab,
      col=col, axes=FALSE, ...)

    if(!is.null(sse)) {
        h <- bounds( which(sse$sse == "H") )
        e <- bounds( which(sse$sse == "E") )

        rect(h[,"start"] , ylim[1], h[,"end"], ylim[1]+4, col="grey20", border=NA)
        rect(e[,"start"] , ylim[1], e[,"end"], ylim[1]+4, col="grey80", border=NA)

        rect(xlim[1], h[,"start"], xlim[1]+4, h[,"end"], col="grey20", border=NA)
        rect(xlim[1], e[,"start"], xlim[1]+4, e[,"end"], col="grey80", border=NA)
    }

    if(!is.null(regs)) {
        r = regs
        rect(r[,"start"] , ylim[2] - 4, r[,"end"], ylim[2], col=col.reg, border=NA)
        rect(xlim[2]-4, r[,"start"], xlim[2], r[,"end"], col=col.reg, border=NA)
    }

    at <- seq(min(x), max(x), length.out=7)
    labs <- resno[ at ]

    axis(1, at = at, labels = labs)
    axis(2, at = at, labels = labs)
    box()




}




