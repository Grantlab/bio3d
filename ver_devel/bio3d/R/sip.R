sip <- function(...)
  UseMethod("sip")

sip.nma <- function(a, b, ...) {
  return(sip.vector(a$fluctuations, b$fluctuations))
}
  

sip.enma <- function(modes, ncore=NULL, ...) {
  if(!inherits(modes, "enma"))
    stop("provide a 'enma' object as obtain from function 'nma.pdbs()'")
  
  ncore <- setup.ncore(ncore, bigmem = FALSE)
  
  if(ncore>1)
    mylapply <- mclapply
  else
    mylapply <- lapply
  
  dims <- dim(modes$fluctuations)
  m <- dims[1]

  mat <- matrix(NA, m, m)
  ##inds <- pairwise(m)
  inds <- rbind(pairwise(m),
                matrix(rep(1:m,each=2), ncol=2, byrow=T))
  
  mylist <- mylapply(1:nrow(inds), function(row) {
    i <- inds[row,1]; j <- inds[row,2];
    out <- list(val=sip.vector(modes$fluctuations[i,], modes$fluctuations[j,]), i=i, j=j)
    return(out)
  })

  for ( i in 1:length(mylist)) {
    tmp <- mylist[[i]]
    mat[tmp$i, tmp$j] <- tmp$val
  }

  mat[ inds[,c(2,1)] ] = mat[ inds ]
  ##diag(mat) <- rep(1, n)
  colnames(mat) <- basename(rownames(modes$fluctuations))
  rownames(mat) <- basename(rownames(modes$fluctuations))
  
  return(round(mat, 6))
}

sip.numeric <- function(...)
  sip.vector(...)

sip.vector <- function(v, w, ...) {
  return(as.numeric(((t(v) %*% w)**2) / ((t(v) %*% v)*(t(w) %*% w))))
}

