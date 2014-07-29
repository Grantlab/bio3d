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

   str <- strsplit(names(grps), "\\.")

   # check for available color code
   ngrps <- rep(NA, length(grps))
   nn = length(unique(grps))
   for(i in 1:nn) {
     inds = which(grps==i)
     nid = as.numeric(sapply(str[inds], "[[", 1))
     cid <- sapply(str[inds], "[[", 2)
     col.ind = lapply(nets[nid], "[[", "community.reindex")
     color = NULL
     for(j in 1:length(col.ind)) {
        if(!is.null(col.ind[[j]])) color = c(color, col.ind[[j]][cid[j]])
     }
     if(length(color) > 1) {
        warning(paste("Ambiguous color in cluster", i))
        color = min(color)
     }
     if(!is.null(color)) ngrps[grps==i] = color 
   }
   tgrp = setdiff(1:nn, ngrps)
   names(tgrp) = unique(grps[is.na(ngrps)])  # new grps
   
   ngrps[is.na(ngrps)] = tgrp[grps[is.na(ngrps)]]

   # Re-coloring for aligned communities
   new.nets <- nets
   max.ncom = max(as.numeric(sapply(str, "[[", 2)))
   ncommunity = rep(list(rep(NA, max.ncom)), length(nets))
   for(i in 1:length(str)) {
     nid <- as.numeric(str[[i]][1])
     cid <- as.numeric(str[[i]][2])
     res <- nets[[nid]]$communities$membership==cid
     V(new.nets[[nid]]$network)$color[res] <- col[ngrps[i]]
     V(new.nets[[nid]]$community.network)$color[cid] <- col[ngrps[i]]
     
     ncommunity[[nid]][cid] <- ngrps[i]
   }
   ncommunity <- lapply(ncommunity, function(x) x[!is.na(x)])

   for(i in 1:length(new.nets)) 
      new.nets[[i]]$community.reindex = ncommunity[[i]] 
   return(list(dismat=dismat, hc=hc, grps=grps, nets=new.nets))
}
