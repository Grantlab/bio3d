`dccm` <-
function(xyz, reference=apply(xyz,2,mean)) {
  
  dotproductmean <- function(xyz, i, j=i) {
    return(  mean(xyz[,atom2xyz(i)] * xyz[,atom2xyz(j)]) )
  }
  
  natm  <- ncol(xyz)/3
  dxyz  <- sweep(xyz, 2, reference)
  dxyz2 <- dxyz^2

  ccmat <- matrix(1, nrow=natm, ncol=natm)
  sqrtdsq <- rep(NA, natm)

  for( i in 1:natm ){
    ##- sqrtdsq[i] = sqrt(dotproductmean(dxyz,i)
    sqrtdsq[i] <- sqrt( mean(dxyz2[,atom2xyz(i)]) )
  }

  ##- Which is quicker? system.time()  
  for(i in 2:natm) {
    for(j in 1:(i-1)) {
      ccmat[i,j] = ccmat[j,i] = dotproductmean(dxyz,i,j)/
        (sqrtdsq[i] * sqrtdsq[j])
    }
  }
 
##  inds <- pairwise(natm)
##  for( i in 1:nrow(inds) ) {
##    ccmat[inds[i,1], inds[i,2]] <- ccmat[inds[i,2], inds[i,1]] <-
##      dotproductmean(dxyz,inds[i,1],inds[i,2])/
##        (sqrtdsq[ inds[i,1] ] * sqrtdsq[ inds[i,2] ])
##  }
  class(ccmat)=c("dccm","matrix")
  return(ccmat)
}

