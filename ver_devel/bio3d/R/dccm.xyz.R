`dccm.xyz` <-
function(x, reference=NULL, grpby=NULL, method=c("pearson", "lmi"), 
         ncore=1, nseg.scale=1, ... ) {
  
  method = match.arg(method)
  
  xyz <- x
  # Parallelized by parallel package (Wed Dec 12 18:36:39 EST 2012)
  ncore <- setup.ncore(ncore)

  if(is.null(reference)) {
     ref = colMeans(xyz)
  } 
  else {
     ref = reference
  }
  dxyz  <- sweep(xyz, 2, ref)

  covmat <- cov(dxyz)
  
  if(!is.null(reference)) {
     # moment instead of covariance
     mxyz <- colMeans(dxyz)
     covmat <- covmat + outer(mxyz, mxyz)
  }

  ccmat <- .cov2dccm(covmat, method = method, ncore = ncore)
 
  if(is.null(grpby)) {
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

# This internal function calculates the N-by-N cross-correlation matrix 
# directly from a 3N-by-3N variance-covariance matrix. 
.cov2dccm <- function(vcov, method = c("pearson", "lmi"), ncore = NULL) {
  method = match.arg(method)
  
  ncore = setup.ncore(ncore)
  if(ncore == 1) mclapply = lapply
  
  x <- vcov
  ccmat = switch(method, 
                 pearson = {
                   n <- nrow(x)
                   np <- pairwise(n/3)
                   
                   d <- sqrt(colSums(matrix(diag(x), nrow=3)))
                   
                   ltv <- mclapply(1:nrow(np), function(i) {
                     i1 <- (np[i, 2] - 1) * 3 + 1
                     i2 <- (np[i, 1] - 1) * 3 + 1
                     sum(diag(x[i1:(i1+2), i2:(i2+2)]))/ # sum of diagnol of submatrix
                       (d[np[i, 2]] * d[np[i, 1]])   # divided by product of standard deviations
                   } )
                   
                   ccmat <- matrix(0, n/3, n/3)
                   ccmat[lower.tri(ccmat)] <- unlist(ltv)
                   
                   # make full matrix
                   ccmat <- ccmat + t(ccmat)
                   diag(ccmat) <- 1
                   ccmat
                 },
                 lmi = {
                   # rm:r-value matrix
                   cm <- x
                   l <- dim(cm)[1]/3
                   rm <- matrix(nrow=l, ncol=l)
                   d <- 3
                   ij <- pairwise(l)
                   
                   # list1: marginal-covariance 
                   list1 <- mclapply(1:l, function(i) det(cm[(3*i-2):(3*i), (3*i-2):(3*i)]) )
                   dm <- unlist(list1)
                   
                   # list2: pair-covariance
                   list2 <- mclapply(1:nrow(ij), function(i) {
                     x <- det(cm[c((3*ij[i,1]-2):(3*ij[i,1]),(3*ij[i,2]-2):(3*ij[i,2])), c((3*ij[i,1]-2):(3*ij[i,1]),(3*ij[i,2]-2):(3*ij[i,2]))])
                     y <- 1/2 * (log(dm[ij[i,1]]) + log(dm[ij[i,2]]) - log(x))
                     (1 - exp(-2 * y / d))^(1/2)
                   }
                   )
                   list2 <- unlist(list2)
                   
                   for (k in 1:nrow(ij)) {
                     rm[ij[k, 1], ij[k, 2]] <- list2[k]
                   }
                   rm[lower.tri(rm)] = t(rm)[lower.tri(rm)]
                   diag(rm) <- 1
                   rm
                 } )
  class(ccmat)=c("dccm","matrix")
  return(ccmat)
}
