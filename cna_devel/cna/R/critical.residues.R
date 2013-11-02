critical.residues <- function(net, method="abs.cij"){

  ## Creates a table with the list of the community-crossing edges
  ## their betweenness value, the nodes forming the edge and their
  ## community membership 

  oops <- require(igraph)
  if (!oops) {
    warning("igraph package missing: Please install, see: ?install.packages")
  }

  if(class(net) != 'cna')
    stop("'net' must be a 'cna' class object")

  ## Function to obtain the cijs of the edges
  edges.cij <- function(cij.matrix,edges){
    values <- NULL
    for(t in 1:dim(edges)[1]){
      values[t] <- round(cij.matrix[as.numeric(edges[t,1]),as.numeric(edges[t,2])],2)
    }
    return(values)
  }

  ## Getting the input
  network.graph <- net$raw.network
  communities <- net$raw.communities

  if(method=="abs.cij"){
    cij.matrix <- net$raw.cij
    diag(cij.matrix) <- 0
  }
  
  ## Edge list generation
  z <- c(1:length(E(network.graph)))

  method.options = c("btwn", "abs.cij")
  method <- match.arg(tolower(method), method.options)
  edges <- matrix(get.edgelist(network.graph), ncol=2)
  
  edge.values <- switch( method,
                       "btwn" = edge.betweenness(network.graph),
                       "abs.cij" =  edges.cij(cij.matrix,edges)
                       )

  edge.list <- cbind(edge.values, edges[,1], edges[,2], z)
  ## create a vector with the crossing edges, their betweenness and the nodes names
  edge.list.crossing.all <- edge.list[crossing(communities,network.graph),] 

  ## Creation of a table with the edge value, the nodes involved and their community membership and the  edge number
  edge.list.crossing.all <- cbind(edge.list.crossing.all[,1:3],communities$membership[edge.list.crossing.all[,2]], communities$membership[edge.list.crossing.all[,3]], edge.list.crossing.all[,4])  

  ## Selection of the edge with the highest btwn/abs.cij per connected community pair
  table.to.screen <- edge.list.crossing.all
  ## Selection of communities
  number.of.comms <- sort(as.numeric(unique(c(table.to.screen[,4], table.to.screen[,5]))))

  highest.value <- NULL
  
  for(y in number.of.comms){
    temporary.table <- NULL
    
    for(t in 1:dim(table.to.screen)[1]){
      if(table.to.screen[t,4] == y){
        temporary.table <- rbind(temporary.table,table.to.screen[t,])
      }
      ## In case the community membership is inverted in the table
      else if(table.to.screen[t,5] == y){  
        invertion <- c(table.to.screen[t,1], table.to.screen[t,3], table.to.screen[t,2], table.to.screen[t,5], table.to.screen[t,4], table.to.screen[t,6])
        names(invertion) <- colnames(table.to.screen)
        temporary.table <- rbind(temporary.table,invertion)
      }
    }

    linked.communities <- sort(as.numeric(unique(temporary.table[,5])))
    for(p in 1:length(linked.communities)){
      same.couple.comms <- NULL

      ## Comparison betweenness values among the same couple of comminities
      for(u in 1:dim(temporary.table)[1]){   
        if(linked.communities[p] == temporary.table[u,5]){
          same.couple.comms <- rbind(same.couple.comms,temporary.table[u,])
        }
      }

      if(method=="btwn"){
        index.highest <- which(same.couple.comms[,2]==max(same.couple.comms[,2]))
      }
      if(method=="abs.cij"){
        index.highest <- which(same.couple.comms[,2]==max(same.couple.comms[,2]))
      }
      highest.crossing <- same.couple.comms[index.highest,]

      highest.value <- rbind(highest.value ,highest.crossing)
    }
  } 

  ## Removal of duplicated lines
  removal.index <- NULL

  for(s in 1:(dim(highest.value)[1]-1)){
    for(z in (s+1):dim(highest.value)[1]){
      if((highest.value[s,4] == highest.value[z,4]) | (highest.value[s,4] == highest.value[z,5])){
        if((highest.value[s,5] == highest.value[z,4]) | (highest.value[s,5] == highest.value[z,5])){
          removal.index <- c(removal.index,s)
        }
      }
    }
  }

  crossing.edges <- highest.value[-removal.index,]

  colnames(crossing.edges) <- c(method, "node1","node2","comm-node1","comm-node2","edge #")
  rownames(crossing.edges) <- NULL
  
  return(crossing.edges)
}
