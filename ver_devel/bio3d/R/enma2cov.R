nma2cov <- function(modes) {
  dims <- dim(modes$U)
  cov <- matrix(0, ncol=dims[1], nrow=dims[1])
  tmpU <- modes$U[,,i]
  tmpL <- modes$L[i,]
  
  for(j in 1:ncol(tmpU) ) {
    cov <- cov + ( (tmpU[,j] %*% t(tmpU[,j])) / tmpL[j])
  }

  return(cov)
}

enma2covs <- function(modes, ncore=4, ...) {
  if(!inherits(modes, "enma"))
    stop("provide a 'enma' object as obtain from function 'nma.pdbs()'")
    
  ncore <- setup.ncore(ncore, bigmem = FALSE)
  
  if(ncore>1)
    mylapply <- mclapply
  else
    mylapply <- lapply
  
  if(!inherits(modes, "enma"))
    stop("provide 'enma' object as obtained from nma.pdbs")
  
  dims <- dim(modes$U.subspace)
  
  mycalc <- function(i, modes) {
    cov <- matrix(0, ncol=dims[1], nrow=dims[1])
    tmpU <- modes$U.subspace[,,i]
    tmpL <- modes$L[i,]

    for(j in 1:ncol(tmpU) ) {
      cov = cov + ( (tmpU[,j] %*% t(tmpU[,j])) / tmpL[j])
    }
    cat(".")
    return(cov)
  }

  covs.list <- mylapply(1:dims[3L], mycalc, modes)
  cat("\n")
  
  covs <- array(0, dim=c(dims[1], dims[1], dims[3]))
  
  for ( i in 1:dims[3L] )
    covs[,,i]=covs.list[[i]]

  ##modes$covs <- covs
  return(covs)
}

.tr <- function(mat) {
  return(sum(diag(mat)))
}

