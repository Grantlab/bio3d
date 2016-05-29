.diag.rtb <- function(H, blocks, xyz, masses, ...) {

  ## mass weighted Hessian
  if(!all(masses == 1)) {
    m <- 1/sqrt(rep(masses, each=3))
    H <- sapply(1:ncol(H), function(i) H[,i] * m[i])
    H <- t( sapply(1:nrow(H), function(i) H[i,] * m[i]) )
  }

  if(nlevels(blocks) == (nrow(H)/3))
    return( eigen(H, ...) )

  MASS <- tapply(masses, blocks, sum)  # Block mass
  xyz <- matrix(xyz, ncol=3, byrow=TRUE)

  ## block center-of-mass
  Rc <- tapply(1:length(blocks), blocks, function(i) {
      colSums(xyz[i,, drop=FALSE ] * masses[i]) / MASS[blocks[i[1]]]
  } )

  ## Gram-Schmidt process for orthonormalization
  schmidt <- function(p) {
    p[, 1] <- p[, 1] / sqrt(sum(p[, 1]^2))
    for(i in 2:ncol(p)) {
      ov <- apply(p[, 1:(i-1), drop=FALSE], 2, function(x)
              crossprod(p[, i], x) * x )
      ov <- rowSums(ov)
      p[, i] <- p[, i] - ov
      p[, i] <- p[, i] / sqrt(sum(p[, i]^2))
    }
    return(p)
  }

  ## RTB projector building (can be parallelized)
  P.blocks <- lapply(1:nlevels(blocks), function(i) {
     myblock <- levels(blocks)[i]
     iatom <- which(blocks %in% myblock)
     natom <- length(iatom)
     m <- masses[iatom]
     M <- MASS[myblock]

     # equilibrium position relative to COM
     r0  <- t( t(xyz[iatom, ]) - Rc[[myblock]] )

     ncol <- ifelse(natom>=3, 6, ifelse(natom>=2, 5, 3))

     P <- matrix(0, nrow=natom*3, ncol=ncol)

     # translation
     P[seq(1, nrow(P), 3), 1] <- sqrt(m/M) #xx
     P[seq(2, nrow(P), 3), 2] <- sqrt(m/M) #yy
     P[seq(3, nrow(P), 3), 3] <- sqrt(m/M) #zz

     if(natom>=3) {

       # rotation
       rr <- rbind(0, -r0[, 3], r0[, 2])
       P[, 4] <- 1/sqrt( sum(m * (r0[, 2]^2 + r0[, 3]^2)) ) *
          ( sqrt(rep(m, each=3)) * as.numeric(rr) )  ## x-axis

       rr <- rbind(r0[, 3], 0, -r0[, 1])
       P[, 5] <- 1/sqrt( sum(m * (r0[, 1]^2 + r0[, 3]^2)) ) *
          ( sqrt(rep(m, each=3)) * as.numeric(rr) )  ## y-axis

       rr <- rbind(-r0[, 2], r0[, 1], 0)
       P[, 6] <- 1/sqrt( sum(m * (r0[, 1]^2 + r0[, 2]^2)) ) *
          ( sqrt(rep(m, each=3)) * as.numeric(rr) )  ## z-axis

       schmidt(P) # Orthonormalization

     } else if(natom==2) {

       # check the line direction and choose the least overlapped axes to rotate.
       ax <- rbind(c(2,3,1), c(3,1,2), c(1,2,3))
       ax <- ax[which.max(abs(r0[1, ])), ]
       rr <- list(rbind(0, -r0[, 3], r0[, 2]),
                  rbind(r0[, 3], 0, -r0[, 1]),
                  rbind(-r0[, 2], r0[, 1], 0))
       # rotation
       P[, 4] <- 1/sqrt( sum(m * (r0[, ax[2]]^2 + r0[, ax[3]]^2)) ) *
          ( sqrt(rep(m, each=3)) * as.numeric(rr[[ax[1]]]) )  ## first-axis

       P[, 5] <- 1/sqrt( sum(m * (r0[, ax[1]]^2 + r0[, ax[3]]^2)) ) *
          ( sqrt(rep(m, each=3)) * as.numeric(rr[[ax[2]]]) )  ## second-axis

       schmidt(P) # Orthonormalization

     } else {

       P

     }
  } )

  ## effective reduced Hessian (can be parallelized)
#  H <- t(P) %*% H %*% P
  H <- lapply(1:nlevels(blocks), function(i) {
    myblock <- unique(blocks)[i]
    iatom <- which(blocks %in% myblock)
    P <- P.blocks[[myblock]]
    crossprod( P, H[atom2xyz(iatom), ] )
  } )
  H <- do.call(rbind, H)

  H <- lapply(1:nlevels(blocks), function(i) {
    myblock <- unique(blocks)[i]
    iatom <- which(blocks %in% myblock)
    P <- P.blocks[[myblock]]
    crossprod( t(H[, atom2xyz(iatom)]), P )
  })
  H <- do.call(cbind, H)

  ei <- eigen(H, ...)

  ## map eigenvector to 3N space
  block.bounds <- c(1, cumsum(sapply(P.blocks[unique(blocks)], ncol))+1)
  vecs <- lapply(1:nlevels(blocks), function(i) {
    myblock <- unique(blocks)[i]
#    iblock <- which(unique(blocks) %in% myblock)
    P <- P.blocks[[myblock]]
    crossprod( t(P), ei$vectors[(block.bounds[i]):(block.bounds[i+1]-1), ] )
  })
  ei$vectors <- do.call(rbind, vecs)

  ei
}


#' @param verbose logical, if TRUE print detailed processing message
#' @inheritParams aanma
#' @rdname aanma
rtb <- function(hessian, pdb, mass=TRUE, nmer=1, verbose=TRUE) {

  H <- hessian

  res <- paste( pdb$atom[, 'chain'],
                pdb$atom[, 'resno'],
                pdb$atom[, 'insert'], sep='_' )

  if(nmer > 1) {
    rl <- rle(res)
    g <- rep(1:ceiling(length(rl$lengths)/nmer), each=nmer)
    tosum <- split(rl$lengths, g[1:length(rl$lengths)])
    rl$lengths <- sapply(tosum, sum)
    rl$values <- rl$values[seq(1, length(g), nmer)]
    res <- inverse.rle(rl)
  }

  blocks <- as.factor(res)

  ## atom mass
  if(isTRUE(mass)) {
      m <- atom2mass(pdb)
  } else if(is.numeric(mass)) {
    m <- mass
  } else {
    m <- rep(1, nrow(pdb$atom))
  }

  if(verbose) {
    cat(" Diagonalizing Hessian with RTB...")
    ptm <- proc.time()
  }

  ei <- .diag.rtb(H, blocks, pdb$xyz, m,  symmetric=TRUE)

  if(verbose) {
    t <- proc.time() - ptm
    cat("\t\tDone in", t[[3]], "seconds.\n")
  }

  return( ei )
}

