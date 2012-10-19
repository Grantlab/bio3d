`dist.xyz` <-
function(a, b=NULL, all.pairs=TRUE){

  ## if 'a' is a vector (or matrix) and 
  ## 'b' is a matrix
  ## compare (each row of) 'a' to all rows in 'b'

  ## if 'a' is a matrix and 'b' is NULL
  ## call 'dist' on 'a'

  ## if 'a' is a vector and 'b' is NULL
  ## make 'a' a 3 col matrix and call 'dist' 


  if(is.vector(a)) {
    a <- matrix(a, ncol=3, byrow=TRUE)
  } else {
    a <- as.matrix(a)
  }
  
  if(is.null(b)) {
    return(as.matrix(dist(a)))
  } else {
    if(is.vector(b)) {
      b <- matrix(b, ncol=3, byrow=TRUE)
    } else {
      b <- as.matrix(b)
    }
  }

  dima <- ncol(a)
  dimb <- ncol(b)
  if(dima != dimb)
    stop("Dimension miss-match of input 'a' and 'b'")
    
  if(dima != 3) {
    warning(paste("input does not have three columns: assuming you want",
                  dima, "dimensional distances")) 
  }

  if(!all.pairs) {
    ## distance between coresponding rows
    d <- rep( NA, max(nrow(a), nrow(b)) )
    ind <- 1:min(nrow(a), nrow(b))
    d[ind] <- sqrt( rowSums((a[ind,] - b[ind,])^2) )
    ##    return( sqrt( rowSums((a - b)^2) ) )
    return(d)
  } else {
    d <- matrix(0, nrow=nrow(a), ncol=nrow(b))
    for(i in 1:nrow(a)){
      d[i,] <- sqrt(colSums((a[i,] - t(b))^2))
    }
    return(d)
  }
}

