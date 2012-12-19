`dccm` <-
function(xyz, reference=apply(xyz,2,mean), ncore=1, nseg.scale=1) {

  # Parallelized by multicore package (Wed Dec 12 18:36:39 EST 2012)
  if(ncore > 1) { 
     require(multicore)
     options(cores = ncore)

     # Issue of serialization problem
     # Maximal number of cells of a double-precision matrix
     # that each core can serialize: (2^31-1-61)/8
     R_NCELL_LIMIT_CORE = 2.68435448e8
     R_NCELL_LIMIT = ncore * R_NCELL_LIMIT_CORE

     if(nseg.scale < 1) {
        warning("nseg.scale should be 1 or a larger integer\n")
        nseg.scale=1
     }
  }
 
  dotproductmean <- function(xyz, i, j=i) {
    return(  mean(xyz[,atom2xyz(i)] * xyz[,atom2xyz(j)]) )
  }
  
  natm  <- ncol(xyz)/3
  dxyz  <- sweep(xyz, 2, reference)
  dxyz2 <- dxyz^2

  ccmat <- matrix(1, nrow=natm, ncol=natm)
  sqrtdsq <- rep(NA, natm)

  for( i in 1:natm ) {
    ##- sqrtdsq[i] = sqrt(dotproductmean(dxyz,i)
    sqrtdsq[i] <- sqrt( mean(dxyz2[,atom2xyz(i)]) )
  }

  if(ncore > 1) {       # Parallelized 
     inds <- pairwise(natm)
     ni <- nrow(inds)
     RLIMIT = R_NCELL_LIMIT
     nDataSeg = floor((ni-1)/RLIMIT)+1
     nDataSeg = floor(nDataSeg * nseg.scale)
     lenSeg = floor(ni/nDataSeg)
     cclist = vector("list", nDataSeg)
     for(i in 1:nDataSeg) {
        istart = (i-1)*lenSeg + 1
        iend = if(i<nDataSeg) i*lenSeg else ni
        cclist[[i]] <- mclapply(istart:iend, function(j)
              dotproductmean(dxyz, inds[j,1], inds[j,2])/
               (sqrtdsq[ inds[j,1] ] * sqrtdsq[ inds[j,2] ]),
           mc.preschedule=TRUE)
     }
     cclist <- unlist(cclist)
     for(i in 1:ni) ccmat[inds[i,1], inds[i,2]] <-
          ccmat[inds[i,2], inds[i,1]] <- cclist[i]
     readChildren()
  } else {       # Single version
     for(i in 2:natm) {
       for(j in 1:(i-1)) {
         ccmat[i,j] = ccmat[j,i] = dotproductmean(dxyz,i,j)/
           (sqrtdsq[i] * sqrtdsq[j])
       }
     }
  }
  ##- Which is quicker? system.time()  
#  for(i in 2:natm) {
#    for(j in 1:(i-1)) {
#      ccmat[i,j] = ccmat[j,i] = dotproductmean(dxyz,i,j)/
#        (sqrtdsq[i] * sqrtdsq[j])
#    }
#  }
 
##  inds <- pairwise(natm)
##  for( i in 1:nrow(inds) ) {
##    ccmat[inds[i,1], inds[i,2]] <- ccmat[inds[i,2], inds[i,1]] <-
##      dotproductmean(dxyz,inds[i,1],inds[i,2])/
##        (sqrtdsq[ inds[i,1] ] * sqrtdsq[ inds[i,2] ])
##  }
  class(ccmat)=c("dccm","matrix")
  return(ccmat)
}

