# Correlation network suboptimal path analysis
#
# Reference 
# Yen, J.Y. (1971) Finding the K Shortest Loopless Paths in a Network.
# Management Science. 17(11):712-716.
cnapath <- function(cna, from, to=NULL, k=10, collapse=TRUE, ncore=NULL, ...) {

  oops <- requireNamespace("igraph", quietly = TRUE)

  if (!oops) 
     stop("igraph package missing: Please install, see: ?install.packages")
  
  if(!inherits(cna, "cna")) {
     if(is.list(cna)) {
        paths <- lapply(cna, cnapath, from=from, to=to, k=k, collapse=collapse, ncore=ncore, ...)
        class(paths) <- c("ecnapath", class(paths))
        return( paths )
     } 
     stop("Input cna is not a 'cna' object")
  }

  ncore = setup.ncore(ncore)

  pairs <- NULL
  if(is.null(to)) {
     if(is.matrix(from)) {
        pairs <- apply(from, 1, function(x) list(from=x[1], to=x[2]))
     } else {
        # all other nodes as sink
        to <- setdiff(seq_along(igraph::V(cna$network)), from)
     }
  } 
  if(is.null(pairs)) {
     from <- unique(from); to <- unique(to)
     for(i in 1:length(from)) 
        for(j in 1:length(to)) pairs <- c(pairs, list(list(from=from[i], to=to[j])))
  }

  ## Initialize progress bar
  pb <- .init.pb(ncore, min=0, max=k*length(pairs))

  ## optimize ncore partition
  ncore.out <- ifelse(ncore < length(pairs), ncore, length(pairs))
  ncore.in <- as.integer(floor(ncore / ncore.out))

  paths <- mclapply(pairs, function(x)
              .cnapath.core(cna=cna, from=x$from, to=x$to, k=k, ncore=ncore.in, 
                 pb=pb, ...), mc.cores=ncore.out)
  
  rm.inds <- sapply(paths, is.null)
  if(sum(rm.inds) == length(paths)) {
     warning(paste("  No path found.\n", 
         "  Please check if the network contains isolated parts!\n\n", sep=""))
     paths <- NULL
  } else {

      if(collapse) {
         paths <- paths[!rm.inds]
         cls <- class(paths[[1]])
         coms <- names(paths[[1]])
         paths <- lapply(coms, function(x) do.call(c, lapply(paths, "[[", x)) )
         names(paths) <- coms 
         class(paths) <- cls
      }
      if(sum(rm.inds)>0) {
         warning(paste("  No path found for some source/sink pairs.\n", 
         "  Please check if the network contains isolated parts!\n\n", sep=""))
      }
  }

  ## Finish progress bar
  .close.pb(pb)

  return(paths)
}

.cnapath.core <- function(cna, from, to, k=1, ncore=1, pb=NULL, ...) {
  
  graph = cna$network

  # which path from the list is the shortest?
  select.shortest.path <- function(variants){
     return( which.min( unlist( lapply( variants, function(x){x$dist} ) ) ) )
  }

  # does a list contain this path?
  contains.path <- function(variants, variant){
     return( any( unlist( lapply( variants, function(x){ isTRUE(all.equal(x$path, variant)) } ) ) ) )
  }

  # first shortest path
  k0 <- suppressWarnings( igraph::get.shortest.paths(graph, from, to, output='both', ...) )
  
  # if no shortest path found, network contains isolated parts.
  if(length(k0$vpath[[1]]) <= 1) {
#     cat("  No path found.\n", 
#         "  Please check if the network contains isolated parts!\n\n", sep="")
     if(!is.null(pb)) .update.pb(pb, step=k)
     return(NULL)
  }

  # number of currently found shortest paths
  kk <- 1
  if(!is.null(pb)) .update.pb(pb)

  # All shortest paths are stored in container A in order
  dist = sum(igraph::E(graph)$weight[k0$epath[[1]]])
  A <- list(list(path=as.integer(k0$vpath[[1]]), epath=as.integer(k0$epath[[1]]), dist=dist))

  # All candidates are stored in container B
  B <- list()

  # until k shortest paths are found
  while(kk < k){
    # take last found shortest path
    last.path <- A[[length(A)]]
    
    tmpB <- mclapply(1:(length(last.path$path)-1), function(i) {    
       spurNode <- last.path$path[i]
       rootPath <- last.path$path[1:i]
       if(i==1) rootePath = NULL
       else rootePath = last.path$epath[1:(i-1)]

       # Remove edges that coincide with the next step from the spur node on
       # those shortest paths stored in A that share the same root path here
       g <- graph
       for(j in 1:length(A)) {
          if(length(A[[j]]$path) > i && isTRUE(all.equal(rootPath, A[[j]]$path[1:i]))) {
             nn = A[[j]]$path[i+1]
             ee = igraph::E(g)[igraph::'%--%'(spurNode, nn)]
             if(length(ee)>0) g <- igraph::delete.edges(g, ee)
          }
       }
       # Remove all edges that link to nodes on the root path (excluding the spur node)
       if(i > 1) {
          for(j in rootPath[-(length(rootPath))]) {
             ee = igraph::incident(g, j) 
             if(length(ee)>0) g <- igraph::delete.edges(g, ee)
          }
       }

       # Suppress warnings because some nodes are intentionally isolated 
       spurPath <- suppressWarnings(igraph::get.shortest.paths(g, spurNode, to, output='both'), ...)

       if(length(spurPath$vpath[[1]]) > 1 ) {
          vpath = c(rootPath, as.integer(spurPath$vpath[[1]][-1]))
          if(!contains.path(B, vpath)) {
             spurPath$epath <- as.integer(igraph::E(graph, path=as.integer(spurPath$vpath[[1]])))
             epath = c(rootePath, spurPath$epath)
             return (list(path=vpath, epath = epath, dist = sum(igraph::E(graph)$weight[epath])) )
          }
       }
       NULL
    }, mc.cores = ncore )
    tmpB <- tmpB[ !sapply(tmpB, is.null) ]
    B <- c(B, tmpB)
    if(length(B) == 0) break
    
    # find shortest candidate
    sp <- select.shortest.path(B)

    # add to A, increase kk, remove shortest path from list of B
    A <- c(A, B[sp])
    kk <- kk + 1
    B <- B[-sp]

    if(!is.null(pb)) .update.pb(pb)
  }

  # stopped before reaching k paths
  if(kk < k) {
    if(!is.null(pb)) .update.pb(pb, k-kk)
    warning("Reaching maximal number of possible paths (", kk, ")")
  }

  out <- list(path=lapply(A, "[[", "path"),  
              epath = lapply(A, "[[", "epath"), 
              dist = sapply(A, "[[", "dist"))
  class(out) <- "cnapath"
  return(out)
}
