#' Align communities from two or more networks
#'
#' Find equivalent communities from two or more networks and re-assign colors
#' to them in a consistent way across networks. A \sQuote{new.membership} vector is
#' also generated for each network, which maps nodes to community IDs that are 
#' renumbered according to the community equivalency.
#'
#' This function facilitates the inspection on the variance of the community
#' partition in a group of similar networks. The original community numbering
#' (and so the colors of communities in the output of \code{plot.cna} and 
#' \code{vmd.cna}) can be inconsistent across networks, i.e. equivalent 
#' communities may display different colors, impeding network comparison. 
#' The function calculates the dissimilarity between all communities and 
#' clusters communities with \sQuote{hclust} funciton. In each cluster, 0 or
#' 1 community per network is included. The color attribute of communities is 
#' then re-assigned according to the clusters through all networks. In addition,
#' a \sQuote{new.membership} vector is generated for each network, which mapps 
#' nodes to new community IDs that are numbered consistently across networks. 
#'
#' @param x,... two or more objects of class \code{cna} (if the numbers of
#'    nodes are different, an alignment \sQuote{fasta} object is required for 
#'    the \code{aln} argument; See below) as obtained from function \code{\link{cna}}. 
#'    Alternatively, a list of \code{cna} objects can be given to \code{x}.
#' @param aln alignment for comparing networks with different numbers of nodes.
#' 
#' @return Returns a list of updated \code{cna} objects.
#'
#' @seealso \code{\link{cna}}, \code{\link{plot.cna}}, \code{\link{vmd.cna}}
#'
#' @examples
#' \donttest{
#'   # Needs MUSCLE installed - testing excluded
#'   if(check.utility("muscle")) {
#'
#'     ## Fetch PDB files and split to chain A only PDB files
#'     ids <- c("1tnd_A", "1tag_A")
#'     files <- get.pdb(ids, split = TRUE, path = tempdir())
#'     
#'     ## Sequence Alignement
#'     pdbs <- pdbaln(files, outfile = tempfile())
#'     
#'     ## Normal mode analysis on aligned data
#'     modes <- nma(pdbs, rm.gaps=TRUE)
#'     
#'     ## Dynamic Cross Correlation Matrix
#'     cijs <- dccm(modes)$all.dccm
#'  
#'     ## Correlation Network
#'     nets <- cna(cijs, cutoff.cij=0.3)
#'
#'     ## Align network communities
#'     nets.aln <- community.aln(nets)
#'
#'     ## plot all-residue and coarse-grained (community) networks
#'     pdb <- pdbs2pdb(pdbs, inds=1, rm.gaps=TRUE)[[1]]
#'     op <- par(no.readonly=TRUE)
#'
#'     # before alignment
#'     par(mar=c(0.1, 0.1, 0.1, 0.1), mfrow=c(2,2))
#'     invisible( lapply(nets, function(x) 
#'        plot(x, layout=layout.cna(x, pdb=pdb, k=3, full=TRUE)[, 1:2], 
#'                full=TRUE)) )
#'     invisible( lapply(nets, function(x) 
#'        plot(x, layout=layout.cna(x, pdb=pdb, k=3)[, 1:2])) )
#'
#'     # after alignment
#'     par(mar=c(0.1, 0.1, 0.1, 0.1), mfrow=c(2,2))
#'     invisible( lapply(nets.aln, function(x) 
#'        plot(x, layout=layout.cna(x, pdb=pdb, k=3, full=TRUE)[, 1:2], 
#'                full=TRUE)) )
#'     invisible( lapply(nets.aln, function(x) 
#'        plot(x, layout=layout.cna(x, pdb=pdb, k=3)[, 1:2])) )
#'
#'     par(op)     
#'   }
#' }
#' @keywords analysis
community.aln <- function(x, ..., aln=NULL) {
    
   ## Check for presence of igraph package
   oops <- requireNamespace("igraph", quietly = TRUE)
   if (!oops) {
      stop("igraph package missing: Please install from CRAN")
   }

   dots <- list(...)
   if(inherits(x, 'cna')) x <- list(x)
   nets <- c(x, dots)
 
   if(length(nets)==1)
     stop('Provide at least two networks.')

   ## Construct dissimilarity matrix
   raw.mat <- lapply(1:length(nets), function(j, aln) {
     net <- nets[[j]]
     ids <- unique(net$communities$membership)
     mat <- NULL
     for(i in ids) {
       r <- t(as.numeric(net$communities$membership==i))
       if(inherits(aln, "fasta")) {
         r2 <- rep(0, ncol(aln$ali))
         r2[!is.gap(aln$ali[j, ])] <- r
         mat <- rbind(mat, r2)
       }
       else {
         mat <- rbind(mat, r)
       }
     }
     rownames(mat) <- ids
     return(mat)
   }, aln)
   
   raw.mat <- do.call(rbind, raw.mat)
   
   ncomms <- sapply(nets, function(net) 
     length(unique(net$communities$membership)))
   r <- rep(1:length(nets), ncomms)
   rownames(raw.mat) <- paste(r, ".", rownames(raw.mat), sep="")
   
   dismat <- dist(raw.mat, method="binary")
   
   ## Clustering
   hc <- hclust(dismat)
   grps <- cutree(hc, h=0.99)
   #plot(hc)

   str <- strsplit(names(grps), "\\.")

   ## check if two or more communities from the same network are in one cluster
   chk <- tapply(str, grps, function(x) sum(duplicated(sapply(x, '[', 1)))>0)
   if(any(unlist(chk))) 
     stop('Two or more communities from the same network are in one cluster.')

   ## rename grps to make sure the group number is assigned based on the rank of
   ## the network in the network list in which the group appears for the first time. 
   ni <- as.numeric(sapply(str, "[[", 1))
   ci <- as.numeric(sapply(str, "[[", 2))
   minc <- tapply(1:length(grps), grps, function(i){
     ind <- order(ni[i], ci[i])[1]
     c(ni[i[ind]], ci[i[ind]])
   }, simplify=FALSE)
   minc <- do.call(rbind, minc)
   inds <- order(minc[, 1], minc[, 2])
   rl <- rle(grps)
   rl$values[] <- match(rl$values, inds)
   grps <- inverse.rle(rl)
   
   
   ## Renumber communities and update 'nets'
   col <- vmd_colors() ## use vmd colors
   for(i in 1:length(nets)) {
     net <- nets[[i]]
     inds <- which(sapply(str, '[[', 1) %in% i)
     mycommid <- as.numeric(sapply(str, '[[', 2))[ inds ]
     mygrps <- grps[inds]
     names(mygrps) <- mycommid
     mygrps <- mygrps[order(mycommid)] # sort new grps according to old commun id

     # update node color
     membership <- as.character(net$communities$membership)
     igraph::V(net$network)$color[] <- as.character(col[mygrps[membership]])

     # generate renumbered membership
     new.membership <- as.numeric(mygrps[membership])
     names(new.membership) <- names(net$communities$membership)
     net$new.membership <- new.membership

     # update community network
     if(!is.null(net$community.network)) {
#       community.cij <- net$community.cij
#       rownames(community.cij) <- mygrps
#       colnames(community.cij) <- mygrps
##       community.cij <- community.cij[order(mygrps), order(mygrps)]
#       community.network <-  igraph::graph.adjacency(community.cij,
#                               mode="undirected",
#                               weighted=TRUE,
#                               diag=FALSE)
       igraph::V(net$community.network)$color <- col[mygrps]
#       igraph::V(community.network)$size <- table(net$communities$membership)
#       net$community.cij <- community.cij
#       net$community.network <- community.network
     }
     nets[[i]] <- net
   }
#   return(list(dismat=dismat, hc=hc, grps=grps, nets=new.nets))

   return( nets )
}
