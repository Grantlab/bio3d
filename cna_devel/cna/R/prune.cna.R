prune.cna <- function(x, edges.min=1, size.min=1) {
  
  ##-- Prune nodes based on number of edges and number of members
  ##     prune.cna(net)
  ##
  
  if(class(x)=="cna") {
    y <- summary.cna(x)
    network=x$clustered.network
  } else {
    warning("Input should be a 'cna' class object as obtained from cna()")
    network=x
    y <- NULL
  }
  ## Identify nodes with less than 'edges.min' to other nodes.
  nodes.inds <- which(degree(network) < edges.min)

  ## Identify nodes with size less than 'size.min'
  ##  cant use V(net$network)$size as these can be scaled
  ##  so we will use the summary information in 'y'
  nodes.inds <- c(nodes.inds, which(y$size < size.min))

  if( length(nodes.inds) == 0 ) {
    cat( "No Nodes Will Removed based on edges.min and size.min values" )
#    output = list(clustered.network=network) #, clustered.communities
    output = x
  } else {
    rm.vs <- V(network)[nodes.inds]
    cat( paste("Removing Nodes:", paste(rm.vs, collapse=", ")),"\n")
    ## Print details of removed with edges
    if(!is.null(y)) {
      w <- cbind(y$tbl[rm.vs,c("id","size")],
                 "edges"=degree(network)[rm.vs],
                 "members"=y$tbl[rm.vs,c("members")])
      write.table(w, row.names=FALSE, col.names=TRUE, quote=FALSE,sep="\t")

      ## Residue raw network
      res2rm <- as.numeric(unlist(y$members[rm.vs]))
      x$raw.communities$membership[res2rm] = NA
    }
  
    d <- delete.vertices(network, rm.vs)
    comms <- edge.betweenness.community(d, directed = FALSE)
    ## Will probably want to keep an edited old community object !!!
    
    output <- list("clustered.network"=d,
                   "clustered.communities"=comms,
                   "raw.network"= x$raw.network,  ## UNCHANGED!!!
                   "raw.communities"=x$raw.communities)
  }

  class(output) = class(x)
  return(output)
}

