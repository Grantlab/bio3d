"com.xyz" <-
  function(xyz, mass=NULL) {
  
  if(is.null(mass))
    mass <- rep(1, times=length(xyz)/3)

  if ((length(xyz)/3) != length(mass))
    stop("com.xyz: unequal lengths")
  
  xyz <- matrix(xyz, ncol=3, byrow=T)
  com <- colSums(xyz * mass) / sum(mass)
  return(com)
}
