"pca.xyz" <-
function(xyz, subset = rep(TRUE, nrow(as.matrix(xyz)))) {
  ## Performs principal components analysis on the given "xyz" numeric data
  ## matrix and return the results as an object of class "pca.xyz"

  xyz <- as.matrix(xyz)
  if (any(!is.finite(xyz)))
    stop("infinite or missing values in x")
  dx <- dim(xyz)
  n <- dx[1]; p <- dx[2]
  if (!n || !p)
    stop("0 extent dimensions")
  
  S    <- var(xyz[subset,])          ## coverance matrix
  mean <- apply(xyz[subset,],2,mean) ## mean structure

  ## eigenvectors ("U") & eigenvalues ("L"): [ U'SU=L ]
  prj  <- eigen(S, symmetric = TRUE)
  L <- prj$values
  U <- prj$vectors

  ## fix negative eigenvalues
  ## (these are very small numbers and should be zero)
  L[L<0]<-0
  sdev <- sqrt(L)

  ## scores of "xyz" on the pc's [ z=U'[x-x.mean] ]
  z <- sweep(xyz,2,mean) %*% (U)

  ## atom-wise loadings (norm of xyz eigenvectors)
  au <- apply(U, 2, function(x) {
    sqrt(colSums(matrix(x^2, nrow=3))) })

  
  class(U)="pca.loadings"

  out <- list(L=L, U=U, z=z, au=au,
              sdev=sdev, mean=mean)

  class(out)="pca"; out
}

