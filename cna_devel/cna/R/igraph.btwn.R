igraph.btwn <-
function(adjmatrix=adjmatrix){

  ## checking igraph uploading
  oops <- require(igraph)
  if (!oops) {
    warning("igraph package missing: Please install, see: ?install.packages")
  }

  ## checking dimentions of network.matrix
  if (dim(adjmatrix)[1] != dim(adjmatrix)[2]) {					
    stop("Input 'adjmatrix' should be a NxN matrix, where N is the number of nodes")
  }
  
  ## creating an igraph network object, edges undirected, edge weights from the matrix
  network <- graph.adjacency(adjmatrix,mode="undirected",weighted=TRUE,diag=FALSE)

  ## edge betweenness calculation with Brandes algorithm
  return(edge.betweenness(network))
}
