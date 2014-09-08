"pca.array" <- function(x) {
  cijs0 <- x
  
  ## Construct the input matrix for PCA
  pcij0 <- t(apply(cijs0, 3, function(x) x[upper.tri(x)]))
  pca.cij0 <- pca.xyz(pcij0, use.svd=T)
  pca.cij0$au <- pca.cij0$U

  return(pca.cij0)
}
