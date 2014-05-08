# align communities of network ensemble and recolored nodes
# number of nodes in the ensemble must be the same
commun.aln <- function(nets, col=vmd.colors()) {
   require(igraph)
   
   #Construct dissimilarity matrix
   raw.mat <- lapply(nets, function(net) {
     ids <- unique(net$communities$membership)
     mat <- NULL
     for(i in ids) {
       r <- t(as.numeric(net$communities$membership==i))
       mat <- rbind(mat, r)
     }
     rownames(mat) <- ids
     return(mat)
   })
   raw.mat <- do.call(rbind, raw.mat)
   ncomms <- unlist(lapply(nets, function(net) 
     length(unique(net$communities$membership))))
   r <- rep(1:length(nets), ncomms)
   rownames(raw.mat) <- paste(r, ".", rownames(raw.mat), sep="")
   
   dismat <- dist(raw.mat, method="binary")
   
   # Clustering
   hc <- hclust(dismat)
   grps <- cutree(hc, h=0.99)
   #plot(hc)

   # Re-coloring for aligned communities
   str <- strsplit(names(grps), "\\.")
   new.nets <- nets
   for(i in 1:length(str)) {
     nid <- as.numeric(str[[i]][1])
     cid <- as.numeric(str[[i]][2])
     res <- nets[[nid]]$communities$membership==cid
     V(new.nets[[nid]]$network)$color[res] <- col[grps[i]]
     V(new.nets[[nid]]$community.network)$color[cid] <- col[grps[i]]
#     new.nets[[nid]]$communities$membership[res] <- grps[i]
   }
   
   return(list(dismat=dismat, hc=hc, grps=grps, nets=new.nets))
}
