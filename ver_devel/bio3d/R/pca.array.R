"pca.array" <- function(x) {
  ## Construct the input matrix for PCA
  pcij0 <- t(apply(x, 3, function(y) y[upper.tri(y)]))
  pca.cij0 <- pca.xyz(pcij0, use.svd=TRUE)
  ##pca.cij0$au <- pca.cij0$U
  pca.cij0$au <- NULL
  return(pca.cij0)
}
