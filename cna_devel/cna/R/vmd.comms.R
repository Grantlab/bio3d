vmd.comms <-
function(vmd.out.file=vmd.out.file,adjmatrix=adjmatrix,vnames=colnames(adjmatrix)){

  oops <- require(igraph)
  if (!oops) {
    warning("igraph package missing: Please install, see: ?install.packages")
  }
 
  if(dim(adjmatrix)[1] != dim(adjmatrix)[2]){
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
  
  contract.matrix <- function(cij.abs, membership,## membership=comms$membership,
                              collapse.method="max"){
    
    ## Function to collapse a NxN matrix to an mxm matrix
    ##  where m is the communities of N. The collapse method
    ##  can be one of the 'collapse.options' below
    collapse.options=c("max", "median", "mean", "trimmed")
    collapse.method <- match.arg(tolower(collapse.method), collapse.options)

    ## Fill a 'collapse.cij' community by community matrix
    node.num <- max(membership)
    collapse.cij <- matrix(0, nrow=node.num, ncol=node.num)
    inds <- pairwise(node.num)

    for(i in 1:nrow(inds)) {
      comms.1.inds <- which(membership==inds[i,1])
      comms.2.inds <- which(membership==inds[i,2])
      submatrix <- cij.abs[comms.1.inds, comms.2.inds]

      ## Use specified "collapse.method" to define community couplings
      collapse.cij[ inds[i,1], inds[i,2] ] = switch(collapse.method,
                  max = max(submatrix),
                  median = median(submatrix),
                  mean = mean(submatrix),
                  trimmed = mean(submatrix, trim = 0.1))
    }

    ## Copy values to lower triangle of matrix and set class and colnames
    collapse.cij[ inds[,c(2,1)] ] = collapse.cij[ inds ]
    colnames(collapse.cij) <- 1:ncol(collapse.cij)
    class(collapse.cij) <- c("dccm", "matrix")
    return(collapse.cij)
  }

  
  out.file <- readLines(vmd.out.file,warn=FALSE)

  ## Creating a vector with residue community membership
  number.comms <- length(grep("The residues in community",out.file)) 
  comms.raw <- out.file[grep("The residues in community",out.file)]
  comm.membership <- NULL

  residues.all <- NULL
  for(comms.index in 1:number.comms){
    residues.raw <- sub(paste("The residues in community ",comms.index," are: ",sep=""), "",comms.raw[comms.index])
    residues <- as.numeric(unlist(strsplit(residues.raw," "))) + 1
 
    residues.all <- c(residues.all, residues)
  }

  node.num <- length(residues.all) 


  comm.membership.raw <- c(1:node.num)
  for(comms.index in 1:number.comms){
    residues.raw <- sub(paste("The residues in community ",comms.index," are: ",sep=""), "",comms.raw[comms.index])
    residues <- as.numeric(unlist(strsplit(residues.raw," "))) + 1

    element <- 1

    for(checking.membership in 1:node.num){
      if((comm.membership.raw[checking.membership] == residues[element]) & (checking.membership <= residues[length(residues)])){
        comm.membership[checking.membership] <- comms.index
        element <- element + 1
      }
    }
    
  }

  ## Creating a table with the critical nodes
  number.crit.nodes <- length(grep("The highest score in edge connectivities between communities",out.file)) 
  crit.raw <- out.file[grep("The highest score in edge connectivities between communities",out.file)]

  btwn <- NULL
  comms <- matrix(NA,ncol=2,nrow=number.crit.nodes)

  for(btwn.index in 1:number.crit.nodes){
    btwn.raw <- sub("The highest score in edge connectivities between communities ", "",crit.raw[btwn.index])
    btwn.raw <- sub("with highest edge connectivity is ", "",btwn.raw)
    btwn.raw <- sub("and ", "",btwn.raw)
    btwn.raw <- unlist(strsplit(btwn.raw," "))
    comms[btwn.index,] <- btwn.raw[1:2]
    #comms[btwn.index,2] <- btwn.raw[2]
    btwn[btwn.index] <- btwn.raw[3] 
  }

  crit.nodes.inds <- NULL
  for(i in 1:length(btwn)){
    check <- grep(btwn[i],out.file)
    crit.nodes.inds[i] <- check[length(check)]
  }

  crit.nodes.temp <- matrix(NA,ncol=3,nrow=number.crit.nodes)
  for(i in 1:length(crit.nodes.inds)){
    crit.nodes.temp[i,] <- unlist(strsplit(out.file[crit.nodes.inds[i]]," "))
  }
    
  crit.nodes <- matrix(NA,ncol=5,nrow=number.crit.nodes)
  for(i in 1:dim(crit.nodes.temp)[1]){
    crit.nodes[i,] <- cbind((as.numeric(crit.nodes.temp[i,3])/2), (as.numeric(crit.nodes.temp[i,1])+1), (as.numeric(crit.nodes.temp[i,2])+1), as.numeric(comms[i,1]), as.numeric(comms[i,2]))
  }
  
  colnames(crit.nodes) <- c("betwn","node1","node2","com-mem-1","com-mem-2")

  membership <- comm.membership
  names(membership) <- vnames

  ## Igraph network generation for plotting VMD-DNA result
  raw.network <- graph.adjacency(adjmatrix,
                             mode="undirected",
                             weighted=TRUE,
                             diag=FALSE)
  
  V(raw.network)$color <- membership+1
  weights <- E(raw.network)$weight * 10

  coarse.cij <- contract.matrix(adjmatrix, membership=membership, collapse.method="max")

  clustered.network <- graph.adjacency(coarse.cij,
                             mode="undirected",
                             weighted=TRUE,
                             diag=FALSE)

  V(clustered.network)$colors <- vmd.colors(length(V(clustered.network)))
  V(clustered.network)$size <- table(membership)

  ## To match summary.cna input names
  raw.communities <- NULL
  raw.communities$membership <- membership

  ## WARNING: THIS IS AN OBJECT JUST TO MATCH summary.cna INPUT NAMES
  ## VMD DOES NOT PERFORM A CLUSTERING OF THE COMMUNITIES
  clustered.communities <- NULL
  clustered.communities$membership <- rep("1",max(membership))
  
  vmd.comms <- list("raw.communities"=raw.communities, "crit.nodes"=crit.nodes,"raw.network"=raw.network,"raw.weights"=weights, "clustered.network"=clustered.network, "clustered.communities"=clustered.communities)

  class(vmd.comms) <- "cna"
  return(vmd.comms)
 

}
