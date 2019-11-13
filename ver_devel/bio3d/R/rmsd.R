"rmsd" <-
function(a, b=NULL,
                 a.inds=NULL,
                 b.inds=NULL,
                 fit=FALSE,
                 ncore=1,
                 nseg.scale=1) {

  # Calculate the RMSD between all rows of 'a' or between
  # the single structure 'a' and the one or more structures
  # contained in 'b'
   
  # Parallelized by parallel package -Wed Dec 12 11:15:20 EST 2012
  # nseg.scale - to resolve the memory problem of using multicore
  ncore <- setup.ncore(ncore)
  if(ncore > 1) {
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

  ## from 'select' object to indices
  if(is.select(a.inds)) a.inds <- a.inds$xyz
  if(is.select(b.inds)) b.inds <- b.inds$xyz

  ## function to fetch xyz and ids from input
  getxyz <- function(x) {
      xyz <- NULL; ids <- NULL;
      
      if(is.pdbs(x)) {
          ids <- basename.pdb(as.character(x$id))
          xyz <- x$xyz
      }
      
      if(is.pdb(x)) {
          xyz <- x$xyz
          if(!is.null(rownames(xyz)))
              ids <- rownames(xyz)
      }
      
      if(is.matrix(x) | is.vector(x)) {
          xyz <- x
          if(!is.null(rownames(xyz)))
              ids <- basename.pdb(as.character(rownames(xyz)))
      }
      
      return(list(xyz=xyz, ids=ids))
  }
    
  ## set xyz and ids for a
  ax <- getxyz(a)
  a=ax$xyz; ids <- ax$ids;

  ## set xyz and ids for a
  if(!is.null(b)) {
      bx <- getxyz(b)
      b=bx$xyz; ids <- bx$ids;
  }

  if( is.null(a.inds) && is.null(b.inds) ) {
    a.inds <- gap.inspect(a)$f.inds

    if(!is.null(b)) { 
      a.inds <- intersect(a.inds, gap.inspect(b)$f.inds)
    }
    b.inds <- a.inds
    if(length(a.inds) != length(a)) {
      warning(paste("No indices provided, using the",
                  length(a.inds)/3,  "non NA positions\n"))
    }
  }

  if (is.null(a.inds)) a.inds <- gap.inspect(a)$f.inds
  if (is.null(b.inds) && !is.null(b)) b.inds <- gap.inspect(b)$f.inds

  if(is.vector(a) || nrow(a)==1) {
    if(is.null(b)) {
      stop("No comparison can be made, input was only a single vector 'a'")
    }
  } else {
    if(is.null(b)) {
      # Pair Wise Matrix 'a'
      if( any(is.na(a[,a.inds])) ) {
        stop("error: NA elements present in selected set")
      }
      nseq=nrow(a)
      inds=pairwise(nseq)
      ni <- nrow(inds)

      if(ncore > 1){          # Parallelized
         RLIMIT = R_NCELL_LIMIT
         nDataSeg = floor((ni-1)/RLIMIT)+1
         nDataSeg = floor(nDataSeg * nseg.scale)
         lenSeg = floor(ni/nDataSeg)
         s = vector("list", nDataSeg)
         for(i in 1:nDataSeg) {
            istart = (i-1)*lenSeg + 1
            iend = if(i<nDataSeg) i*lenSeg else ni 
            s[[i]] <- mclapply(istart:iend, function(j) {
                 x <- a[inds[j,1],a.inds]
                 y <- a[inds[j,2],a.inds]
                 if(fit) {
                   y <- fit.xyz(fixed=x, mobile=y, fixed.inds=1:length(x))
                 }
                 sqrt(sum((x-y)^2)/(length(a.inds)/3))
               },
               mc.preschedule=TRUE)
         }
         s <- unlist(s)
      } else {               # Single version
         s <- rep(NA, ni)
         for(i in 1:ni) {
           x <- a[inds[i,1],a.inds]
           y <- a[inds[i,2],a.inds]
           if(fit) {
             y <- fit.xyz(fixed=x, mobile=y, fixed.inds=1:length(x))
           }
           s[i]<-sqrt(sum((x-y)^2)/(length(a.inds)/3))
         }
      }

      sm <- matrix(0, ncol=nseq,nrow=nseq)
      sm[inds]<-s
      if(nseq==2) {
        sm[inds[,2], inds[,1]]<-s
      } else {
        sm[inds[,c(2,1)]]<-s
      }

      if(!is.null(ids)) {
          rownames(sm) <- ids
          colnames(sm) <- ids
      }
      return(round(sm,3))

    } else {
      stop("input 'b' supplied:
              Therefore 'a' should be a single vector for comparison with 'b'.
               Not a matrix or list")
    }
  }

  if(is.vector(b) || nrow(b)==1) {
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
    if(is.matrix(b) && (nrow(b) > 1)) {
      if (length(a.inds) != length(b.inds)) {
          stop("dimension mismatch:
                  a[a.inds] and b[,b.inds] should be the same length")
      }
      if( any(is.na(a[a.inds])) ||
         any(is.na(b[,b.inds])) ) {
        stop("error: NA elements present in selected set")
      }
      if(fit) {
         # Parallelized / single version
         b <- fit.xyz(fixed=a, mobile=b,
                     fixed.inds=a.inds, mobile.inds=b.inds, 
                     ncore=ncore, nseg.scale=nseg.scale)
      }

      if(ncore > 1){            # Parallelized
         RLIMIT = R_NCELL_LIMIT
         nDataSeg = floor((nrow(b)-1)/RLIMIT)+1
         nDataSeg = floor(nDataSeg * nseg.scale)
         lenSeg = floor(nrow(b)/nDataSeg)
         irmsd = vector("list", nDataSeg)
         for(i in 1:nDataSeg) {
            istart = (i-1)*lenSeg + 1
            iend = if(i<nDataSeg) i*lenSeg else nrow(b)
            irmsd[[i]] <- mclapply(istart:iend, function(j) 
                  sqrt( sum((b[j,b.inds]-a[a.inds])^2)/(length(a[a.inds])/3) ),
               mc.preschedule=TRUE)
         }
         irmsd <- unlist(irmsd)
         if(!is.null(ids)) names(irmsd) <- ids
      } else {                  # Single version
         irmsd <- sqrt( apply((apply(b[,b.inds],1,"-",a[a.inds])^2),2,sum)/(length(a[a.inds])/3) )
         if(!is.null(ids)) names(irmsd) <- ids
      }
      return(round(irmsd,3))
    }
  }
}
