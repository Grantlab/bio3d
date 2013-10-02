"seq.identity" <-
function( alignment , normalize=TRUE, ncore=1, nseg.scale=1) {

  # Parallelized by multicore package (Sun Jul  7 17:35:38 EDT 2013)
  if(ncore > 1) {
     oops <- require(multicore)
     if(!oops)
        stop("Please install the multicore package from CRAN")

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
 
  if(is.list(alignment)) alignment <- alignment$ali
  alignment[is.gap(alignment)] = NA

  ide <- function(x, y) {
#### Edit by Heiko Strathmann
#### Wed Aug  4 10:48:16 PDT 2010
#### Fix for bug with all gap sequences
    r <- sum(x==y, na.rm=TRUE)
    t <- sum(complete.cases(cbind(x,y)))
    if (normalize && t != 0) {
      r <- r/t
    }
##################################
    return( round(r, 3) )
  }

  nseq <- nrow(alignment)
  inds <- pairwise( nseq )
  ni <- nrow(inds)

  if(ncore > 1) {
     RLIMIT = R_NCELL_LIMIT
     nDataSeg = floor((ni-1)/RLIMIT) + 1
     nDataSeg = floor(nDataSeg * nseg.scale)
     lenSeg = floor(ni/nDataSeg)
     s = NULL
     for(i in 1:nDataSeg) {
        istart = (i-1)*lenSeg + 1
        iend = if(i<nDataSeg) i*lenSeg else ni
        s <- c(s, mclapply(istart:iend, function(j) {
           ide(alignment[inds[j,1],], alignment[inds[j,2],])
          }) )
     }
     s <- unlist(s)
     readChildren()
  } else {
     s <- rep(NA, ni)
   
     for(i in 1:ni) {
       s[i]<-ide(alignment[inds[i,1],], alignment[inds[i,2],])
     }
  }
  ## make 's' into matrix 'sm'
  sm <- matrix(1, ncol=nseq,nrow=nseq)
  sm[inds]<-s
  if(nseq==2) {
    sm[inds[,2], inds[,1]]<-s
  } else {
    sm[inds[,c(2,1)]]<-s 
  }
  return(sm) # ide matrix
}

