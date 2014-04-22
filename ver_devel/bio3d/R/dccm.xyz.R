`dccm.xyz` <-
function(x, reference=NULL, grpby=NULL, ncore=1, nseg.scale=1, ... ) {
  xyz <- x
  # Parallelized by multicore package (Wed Dec 12 18:36:39 EST 2012)
  ncore <- setup.ncore(ncore)

  if(is.null(reference)) ref = colMeans(xyz)
  else ref = reference
  dxyz  <- sweep(xyz, 2, ref)

  covmat <- cov(dxyz)
  
  if(!is.null(reference)) {
     # moment instead of covariance
     mxyz <- colMeans(dxyz)
     covmat <- covmat + outer(mxyz, mxyz)
  }
  n <- nrow(covmat)
  np <- pairwise(n/3)

  if(ncore > 1) 
     mylapply <- mclapply
  else
     mylapply <- lapply

  ltv <- mylapply(1:nrow(np), function(x) {
     i1 <- (np[x, 2] - 1) * 3 + 1
     i2 <- (np[x, 1] - 1) * 3 + 1
     sum(diag(covmat[i1:(i1+2), i2:(i2+2)]))
  } )
#ltv=list(rep(1, nrow(np)))
  ccmat <- matrix(0, n/3, n/3)
  ccmat[lower.tri(ccmat)] <- unlist(ltv)
  ccmat <- ccmat + t(ccmat)
  diag(ccmat) <- colSums(matrix(diag(covmat), nrow=3))
  d <- sqrt(diag(ccmat))
  d <- outer(d, d)
  ccmat <- ccmat / d

  if(is.null(grpby)) {
    class(ccmat)=c("dccm","matrix")
    return(ccmat)
  } else {
    ##- Group by concetive numbers in 'grpby'
    if( ncol(xyz) != (length(grpby)*3) )
      stop("dimension miss-match in 'xyz' and 'grpby', check lengths")
    
    ##- Bounds of 'grpby' numbers
    inds <- bounds(grpby, dup.inds=TRUE)
    nres <- nrow(inds)

    ##- Per-residue matrix
    m <- matrix(, ncol=nres, nrow=nres)
    ij <- pairwise(nres)

    ##- Max (absolute value) per residue
    for(k in 1 : nrow(ij) ) {
      m[ij[k,1],ij[k,2]] <-
        min( ccmat[ (inds[ij[k,1],"start"]:inds[ij[k,1],"end"]),
                  (inds[ij[k,2],"start"]:inds[ij[k,2],"end"])],
            na.rm=TRUE )
      tmax <- max( ccmat[ (inds[ij[k,1],"start"]:inds[ij[k,1],"end"]),
                  (inds[ij[k,2],"start"]:inds[ij[k,2],"end"])],
            na.rm=TRUE )
      if(tmax > abs(m[ij[k,1],ij[k,2]])) m[ij[k,1],ij[k,2]] = tmax 
    }
#    if( !mask.lower )
    m[lower.tri(m)] = t(m)[lower.tri(m)]
    diag(m) <- 1

    class(m)=c("dccm","matrix")
    return(m)
  }
}

