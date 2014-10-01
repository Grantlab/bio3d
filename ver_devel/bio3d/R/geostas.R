geostas <- function(xyz, amsm=NULL, k=3, method="pairwise", fit=TRUE, ...) {

#  if(all(c(!is.matrix(xyz), !inherits(xyz, "pdbs"), !is.pdb(xyz))))
  if(! (is.matrix(xyz) || inherits(xyz, "pdbs") || is.pdb(xyz)) )
    stop(paste("'xyz' should be a trajectory matrix, a 'pdb' object, or\n\t",
               "a list object as obtained from 'read.fasta.pdb'"))
  
  if(!method %in% c("pairwise", "columnwise"))
    stop("'method' should be 'pairwise' or 'columnwise'")
  
  if(inherits(xyz, "pdbs")) {
    pdbs <- xyz
    gaps.pos <- gap.inspect(pdbs$xyz)
    xyz <- pdbs$xyz[, gaps.pos$f.inds]
  }

  if(is.pdb(xyz)) {
    if(!is.matrix(xyz$xyz))
      stop("incompatible input. provide a multimodel 'pdb' object")
    
    pdb <- xyz
    ca.inds <- atom.select(pdb, 'calpha')
    xyz <- pdb$xyz[,ca.inds$xyz]
  }

  ##if(fit) {
  ##  inds <- seq(1, ncol(xyz))
  ##  xyz <- fit.xyz(xyz[1,], xyz, fixed.inds=inds, mobile.inds=inds)
  ##}
  
  if(fit && is.null(amsm)) {
    core <- core.find(xyz)
    fit.inds <- core$c1A.xyz

    xyz <- fit.xyz(xyz[1,], xyz, fixed.inds=fit.inds, mobile.inds=fit.inds)
    if(is.null(fit.inds))
      warning("core indices not found. fitting to all atoms")
  }
  else {
    fit.inds <- NULL
  }
  
  if(is.null(amsm)) {
    amsm <- amsm.xyz(xyz, ...)
  }
  else {
    if(!all(dim(amsm)==ncol(xyz)/3))
      stop("dimension mismatch ('xyz' and 'amsm')")
  }

  if(method=="pairwise") {
    cm <- 1-amsm
    diag(cm) <- 0
  }
  else {
    cm.tmp  <- normalize.vector(amsm)
    cm <- 1 - apply(cm.tmp, 2, function(x,y) x %*% y, cm.tmp)
  }

  
  ## H-clust
  #dis    <- dist(cm)
  #hc     <- hclust(dis)
  #grps   <- cutree(hc, k = k) 

  
  ## H-clust method: 'ward'
  #d    <- dist(cm, method = "euclidean") # distance matrix
  #fit  <- hclust(d, method="ward") 
  #grps <- cutree(fit, k=k) # cut tree into 5 clusters


  ## k-means clustering
  grps <- kmeans(cm, k)$cluster

  ## return also indices for the identified domains?
  out <- list(grps=grps, amsm=amsm, fit.inds=fit.inds)
  return(out)
}
