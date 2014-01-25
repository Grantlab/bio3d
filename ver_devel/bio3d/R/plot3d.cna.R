plot3d.cna <- function(x,
                       pdb = NULL,
                       node.size = NULL, ##
                       weights=NULL,
                       layout = layout.pdb(pdb, x),
                       col = NULL,
                       ...){
  ##
  ## Plot a cna network graph in 3D with RGL
  ##    plot3d.cna(net, pdb)
  ##  Or just
  ##    rglplot(net$clustered.network)
  ##
  
  ## Check if x has 'cna' class
  if(!"cna" %in% class(x)){
    stop("Input 'x' object must be a 'cna' class object")
  }

  oops <- require(igraph)
  if (!oops) {
    warning("igraph package missing: Please install, see: ?install.packages")
  }

  if(is.null(weights)){
    weights <- E(x$community.network)$weight
    
    if(is.null(x$call$minus.log)){
      weights <- exp(-weights)
    }
    else{
      if(x$call$minus.log){
        weights <- exp(-weights)
      }
    }
     weights <- weights*10
  }
  
  ## Obtain the plot coords...
  if(!is.null(pdb) && is.null(layout)) {
    cat("Obtaning layout from PDB structure\n")
    layout = layout.pdb(pdb, x)
  }
  if(is.null(pdb) && is.null(layout)) {
    cat("Obtaning guestimated layout with fruchterman.reingold\n")
    layout <- layout.fruchterman.reingold(x$community.network, weights=weights)
  }
  if(dim(layout)[2] != 3){
    stop("Input 'layout' must be an Nx3 matrix, where N is the number of communities")
  }
  
  rglplot(x$community.network,
          edge.width = weights,
          layout = layout,
          vertex.size = node.size,
          vertex.color <- col)
}

