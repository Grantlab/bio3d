
covsoverlap <- function(...)
  UseMethod("covsoverlap")

covsoverlap.enma <- function(enma, ncore=4, subset=NULL, ...) {
  if(!inherits(enma, "enma"))
    stop("provide a 'enma' object as obtain from function 'nma.pdbs()'")

    ncore <- setup.ncore(ncore, bigmem = FALSE)

  if(ncore>1)
    mylapply <- mclapply
  else
    mylapply <- lapply

  cat("Calculating pairwise covariance overlap coefs")

  m <- dim(enma$U.subspace)[3]
  mat <- matrix(NA, m, m)
  ##inds <- pairwise(m)
  inds <- rbind(pairwise(m),
                matrix(rep(1:m,each=2), ncol=2, byrow=T))
  
  mylist <- mylapply(1:nrow(inds), function(row) {
    i <- inds[row,1]; j <- inds[row,2];
    val <- covsoverlap.matrix(enma$U.subspace[,,i],
                              enma$U.subspace[,,j],
                              enma$L[i, ],
                              enma$L[j, ],
                              subset=subset)
               
    out <- list(val=val, i=i, j=j)
    cat(".")
    return(out)
  })
  
  for ( i in 1:length(mylist)) {
    tmp <- mylist[[i]]
    mat[tmp$i, tmp$j] <- tmp$val
  }

  mat[ inds[,c(2,1)] ] = mat[ inds ]
  ##diag(mat) <- rep(1, n)

  rownames(mat) <- basename(rownames(enma$fluctuations))
  colnames(mat) <- basename(rownames(enma$fluctuations))
  
  cat("\n")
  return(round(mat, 6))
  
}



covsoverlap.matrix <- function(Ua, Ub, La, Lb, subset=NULL) {
  if(any(missing(Ua), missing(Ub), missing(La), missing(Lb)))
    stop("provide eigenvectors and eigenvalues")
  
  if(!is.null(subset)) {
    if(subset>ncol(Ua))
      subset <- ncol(Ua)
    
    Ua <- Ua[,1:subset]
    Ub <- Ub[,1:subset]
    La <- La[1:subset]
    Lb <- Lb[1:subset]
  }

  sumb <- 0
  for( k in 1:ncol(Ua) ) {
    
    tmp <- sqrt(La[k] * Lb)
    overlap <- c((t(Ua[,k]) %*% Ub)**2)
    sumb <- sumb + sum( tmp * overlap )
  }
  
  return(1 - ( sum(La + Lb) - 2 *sumb ) / sum(La + Lb))
}

