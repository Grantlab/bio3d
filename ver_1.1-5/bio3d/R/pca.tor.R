"pca.tor" <-
function(data, subset = rep(TRUE, nrow(as.matrix(data))) ) {

  data[data < 0]    <- data[data < 0] + 360
  data[data < -180] <- data[data < -180] + 360
  data[data >  180] <- data[data >  180] - 360
  ##cat("Rescaled (-180 to 180) and corrected for periodicity\n")
  data <- wrap.tor(data)
  return( pca.xyz(data, subset) )
}

