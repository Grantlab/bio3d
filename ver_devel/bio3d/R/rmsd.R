"rmsd" <-
function(a, b=NULL,
                 a.inds=NULL,
                 b.inds=NULL,
                 fit=FALSE) {

  # Calculate the RMSD between all rows of 'a' or between
  # the single structure 'a' and the one or more structures
  # contained in 'b'

  if(is.list(a)) a=a$xyz

  if(is.vector(a)) {
    if(is.null(b)) {
      stop("No comparison can be made, input was only a single vector 'a'")
    }
    lseq   <- length(a)
    if(is.null(a.inds)) a.inds <- 1:lseq
  } else {
    lseq   <- ncol(a)
    if(is.null(a.inds)) a.inds <- 1:lseq
    if(is.null(b)) {
      # Pair Wise Matrix 'a'
      if( any(is.na(a[,a.inds])) ) {
        stop("error: NA elements present in selected set")
      }
      nseq=nrow(a)
      inds=pairwise(nseq)
      ni <- nrow(inds)
      s <- rep(NA, ni)

      for(i in 1:ni) {
        x <- a[inds[i,1],a.inds]
        y <- a[inds[i,2],a.inds]
        if(fit) {
          y <- fit.xyz(fixed=x, mobile=y)
        }
        s[i]<-sqrt(sum((x-y)^2)/(length(a.inds)/3))
      }

      sm <- matrix(0, ncol=nseq,nrow=nseq)
      sm[inds]<-s
      if(nseq==2) {
        sm[inds[,2], inds[,1]]<-s
      } else {
        sm[inds[,c(2,1)]]<-s
      }
      return(round(sm,3))

    } else {
      stop("input 'b' supplied:
              Therefore 'a' should be a single vector for comparison with 'b'.
               Not a matrix or list")
    }
  }

  if(is.vector(b)) {
    if(is.null(b.inds)) b.inds=1:length(b)
    if (length(a.inds) != length(b.inds)) {
      stop("dimension mismatch:
              a[a.inds] and b[b.inds] should be the same length")
    }
    if( any(is.na(a[a.inds])) ||
       any(is.na(b[b.inds])) ) {
      stop("error: NA elements present in selected set")
    }
    if(fit) {
      b <- fit.xyz(fixed=a,mobile=b,fixed.inds=a.inds,mobile.inds=b.inds)
    }

    irmsd <- sqrt(sum((a[a.inds] - b[b.inds])^2)/(length(a[a.inds])/3) )
    return(round(irmsd,3))

  } else {
    if(is.list(b)) {
      b=b$xyz
    }
    if(is.matrix(b)) {
      if(is.null(b.inds)) b.inds=1:ncol(b)
        if (length(a.inds) != length(b.inds)) {
          stop("dimension mismatch:
                  a[a.inds] and b[,b.inds] should be the same length")
        }
      if( any(is.na(a[a.inds])) ||
         any(is.na(b[,b.inds])) ) {
        stop("error: NA elements present in selected set")
      }
      if(fit) {
        b <- fit.xyz(fixed=a, mobile=b,
                     fixed.inds=a.inds,mobile.inds=b.inds)
      }

      irmsd <- sqrt( apply((apply(b[,b.inds],1,"-",a[a.inds])^2),2,sum)/(length(a[a.inds])/3) )
      return(round(irmsd,3))

    }
  }
}
