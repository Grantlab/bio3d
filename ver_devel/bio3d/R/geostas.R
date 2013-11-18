geostas <- function(xyz, amsm=NULL, k=3, ...) {
  ## TODO: include option to input pdbs 3dalign object
  
  if(!is.matrix(xyz))
    stop("'xyz' must be a trajectory matrix as obtained e.g. by\n")
  
  if(is.null(amsm))
    amsm <- amsm.xyz(xyz, ...)

  cm.tmp  <- normalize.vector(amsm)
  cm.norm <- apply(cm.tmp, 2, function(x,y) x %*% y, cm.tmp)

  dis    <- dist(cm.norm)
  hc     <- hclust(dis)
  grps   <- cutree(hc, k = k) 
  
  ##hc1 <- reorder.hclust(clusts, dis)
  ##grps <- cutree(hc1, k = k) 

  ## return also indices for the identified domains?
  out <- list(grps=grps, amsm=amsm)
  
  return(out)
}
