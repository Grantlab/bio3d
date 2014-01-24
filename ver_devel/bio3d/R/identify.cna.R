identify.cna <- function(x, labels=NULL, cna=NULL, labelshorten=TRUE, ...){

  ## Be carefull with input argument order
  ##  - 'labels' can take any input and screw up priniting
  ##      e.g. if you pass cna as the second argument!
  ## Should this perhaps be able to take just a cna object as input
  ##   - Possible if cna object has layout defined
  ##   - Could take extra layout option for ustom graphs
  ## xy <- plot.cna(net)
  ## d <- identify.cna(xy, labels=net$raw.communities$membership)
  ## d <- identify(xy, cnet=net)
  ## d <- identify(xy, labels=summary(net)$members)
  ##
  ## ToDo: Fix c(87:87) for cluster 8 in summary.cna() and print out etc.
  
  oops <- require(igraph)
  if (!oops) {
    stop("igraph package missing: Please install, see: ?install.packages")
  }

  if(dim(x)[2] != 2){
    stop("'x' object must be a Nx2 numeric matrix")
  }
  
  x.norm <- layout.norm(x, -1, 1, -1, 1)

  if(is.null(cna) && is.null(labels)) {
    inds <- identify(x.norm[,1], x.norm[,2], ...)
    return(inds)
  }
  if(!is.null(cna) && is.null(labels)) {
    ## take labels from cna object
    ## labels <- summary.cna(cna)$tbl$members
    ## labels <- summary.cna(cna)$members
    lab.all <- summary.cna(cna)
    if(labelshorten) {
      labels <- lab.all$members
      names(labels) <- lab.all$id
    } else {
      labels <- lab.all$members
    }
  }

  inds <- identify(x.norm[,1], x.norm[,2], labels, ...)
  return( labels[inds] )
}
