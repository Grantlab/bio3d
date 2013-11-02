vmd.btwn <-
function(btwn.matrix=btwn.matrix){

  if (dim(btwn.matrix)[1] != dim(btwn.matrix)[2]) {					## checking dimentions of btwn.matrix
    stop("Input 'btwn.matrix' should be a NxN matrix, where N is the number of nodes")
  }

  btwn.matrix[upper.tri(btwn.matrix)] <- 0						## Setting half betweenness matrix to 0
  bet.vmd <- btwn.matrix[btwn.matrix>0] / 2						## Divinding the betweenness edge as VMD-plugin counts the paths from node i to j and from node j to i
  return(bet.vmd)
}
