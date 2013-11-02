comms.significance <- function(adjmatrix, comms=comms, vnames=colnames(adjmatrix)) {
  ## Check for igraph
  oops <- require(igraph)
  if (!oops) {
    warning("igraph package missing: Please install, see: ?install.packages")
  }
  if (dim(adjmatrix)[1] != dim(adjmatrix)[2]) {
    stop("Input 'network.matrix' should be a NxN matrix, where N is the number of nodes")
  }

  ## Check vnames/colnames if present
  if( is.null( vnames ) ) {
    vnames <- paste0("res", 1:ncol(adjmatrix))
    colnames(adjmatrix) <- vnames
  } 
  if( length(vnames) != ncol(adjmatrix) ) {
    stop("Length of input vnames and number of cols in adjmatrix do not match") 
  }

  oops <- attributes(comms)
  if(oops$class != "community"){
    stop("comms argument must be an igraph community object") 
  } 
  if(length(comms$membership) != dim(adjmatrix)[1]){
     stop("adjmatrix and comms argument do not match") 
  } 
  ## Make an igraph network object
  network <- graph.adjacency(adjmatrix,
                             mode="undirected",
                             weighted=TRUE,
                             diag=FALSE)

  significance <- NULL

  for(i in 1:max(comms$membership)){
    subset <- V(network)[comms$membership==i]
    subgraph <- induced.subgraph(network, subset)
    in.degrees <- degree(subgraph)
    out.degrees <- degree(network, subset) - in.degrees
    output.wilcox <- wilcox.test(in.degrees, out.degrees, exact=FALSE, paired=TRUE)
      
    significance[i] <- output.wilcox$p.value
  }

  output <- cbind(c(1:max(comms$membership)), significance)
  colnames(output) <- c("Comms", "P.values")
      
  return(output)
}
