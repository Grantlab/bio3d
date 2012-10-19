"seqbind" <-
function(a, b, blank="-") {

  if(is.null(a) | is.null(b)) {
    raw <- rbind(a,b)
    return(raw)
  } else {

    if(is.vector(a)) {
      ca <- length(a); ra <- 1; ma <- FALSE
    } else {
      ca <- ncol(a); ra <- nrow(a); ma <- TRUE
    }
    if(is.vector(b)) {
      cb=length(b); rb <- 1; mb <- FALSE
    } else {
      cb <- ncol(b); rb <- nrow(b); mb <- TRUE
    }
    
    if(ca!=cb) {
      c <- which.min(c(ca,cb)) 
      if(c==1) {
        if(ma) {
          a <- cbind( a,matrix(blank,ncol=(cb-ca),nrow=nrow(a)) )
        } else {
          a <- c( a, rep(blank,(cb-ca)) )
        }
      } else {
        if(mb) {
          b <- cbind( b,matrix(blank,ncol=(ca-cb),nrow=nrow(b)) )
        } else {
          b <- c( b, rep(blank,(ca-cb)) )
        }
      }
    }
    raw <- rbind(a,b)
    return(raw)
  }
}

