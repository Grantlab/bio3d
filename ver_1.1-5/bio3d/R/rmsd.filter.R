"rmsd.filter" <-
function(xyz=NULL, rmsd.mat=NULL, cutoff=0.5, fit=TRUE, verbose=TRUE) {

  # k<-rmsd.filter(xyz=pdbs$xyz, cutoff=0.5)
  # k<-rmsd.filter(rmsd.mat=k$rmsd.mat, cutoff=2.0)
  
  if(is.null(rmsd.mat)) {
    if(is.null(xyz))
      stop("Must provide either an alignment 'aln' or identity matrix 'ide'")
    if(is.list(xyz))
      xyz=xyz$xyz
    gaps<-gap.inspect(xyz)
    rmsd.mat <- rmsd( xyz[,gaps$f.inds],fit=fit )
  }
  
  r.d  <- as.dist(rmsd.mat)
  tree <- hclust(r.d)

  h = cutoff
  n <- nrow(tree$merge) + 1
  k <- integer(length(h))
  k <- n + 1 - apply(outer(c(tree$height, Inf), h, ">"),2, which.max)
  if(verbose)
    cat("ide.filter(): N clusters @ cutoff = ", k, "\n")
  
  ans <- as.vector(.Call("R_cutree", tree$merge, k, PACKAGE = "stats"))
  
  cluster.rep <- NULL
  for(i in 1:k) {
    ind <- which(ans==i)
    if (length(ind) == 1) {
      cluster.rep <- c(cluster.rep, ind)
    } else {
      cluster.rep <- c(cluster.rep,
                     ind[ which.min( colSums(rmsd.mat[ind,ind]) ) ])
    }
  }
  return(list(ind=cluster.rep, tree=tree, rmsd.mat=rmsd.mat))
}

