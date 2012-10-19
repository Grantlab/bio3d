"identity" <-
function( alignment , normalize=TRUE) {
  
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
  s <- rep(NA, ni)
  
  for(i in 1:ni) {
    s[i]<-ide(alignment[inds[i,1],], alignment[inds[i,2],])
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

