vec2color <- function(vec, pal=c("blue", "green", "red"), n=30) {
  ##-- Define a color scale from a numeric vector
  require(classInt)
  return( findColours(classIntervals(vec, n=n, style="equal"), pal) )
}

normalize.cij <- function(cij, factor = NULL, mag = 2, cutoff = 0) {
   scale <- mag
   if("all.dccm" %in% names(cij)) cij <- cij$all.dccm
   if(is.list(cij))
      cij <- array(unlist(cij), dim = c(dim(cij[[1]]), length(cij)))
   if(!is.array(cij)) stop("Input format should be 2D/3D array, ensmb, or list")
   if(is.na(dim(cij)[3L]))
      cij <- array(cij, dim = c(dim(cij), 1))
   if(is.null(factor))
      factor = factor(1:dim(cij)[3L])
   factor <- as.factor(factor)

   cij <- tapply(1:dim(cij)[3L], factor, function(i)
       rowMeans(cij[,,i, drop=FALSE], dims=2) )
   pcij <- do.call(rbind, lapply(cij, function(x) x[upper.tri(x)]))
   
   # Variance
   vars <- apply(pcij, 2, var)
   diff = max(vars, na.rm=TRUE) - min(vars, na.rm=TRUE)
   if(diff > 0) 
      colors <- (vars - min(vars, na.rm=TRUE))/ diff * scale
   else
      colors <- rep(0, length(vars))
 
   colors <- matrix(rep(colors, nrow(pcij)), nrow=nrow(pcij), byrow=TRUE)
   colors[pcij <= cutoff] <- NA

   # Weight
#   tmax <- apply(pcij, 2, max)
   pcij[pcij <= cutoff] <- NA
   vmin <- min(pcij, na.rm = TRUE)
   vmax <- max(pcij, na.rm = TRUE)
   pcij <- (pcij - vmin) / (vmax - vmin) * 0.9 + 0.05

   pcij[is.na(pcij)] <- 0

   ncij <- apply(pcij, 1, function(x) {
      xx = cij[[1]]
      xx[upper.tri(xx)] <- x
      xx[lower.tri(xx)] <- t(xx)[lower.tri(xx)]
      return(list(xx))
   } )

   ncolor <- apply(colors, 1, function(x) { 
      xx = cij[[1]]
      xx[upper.tri(xx)] <- x
      xx[lower.tri(xx)] <- t(xx)[lower.tri(xx)]
      return(list(xx))
   } )

   cij <- do.call("c", ncij)
   colors <- do.call("c", ncolor)
   if(length(cij) == 1) colors = NULL
   return( list(cij=cij, color=colors) )
}

# None: minus.log of max cij (default of cna()); if member=NULL and col=NULL, no change
# Mean: take mean cij of top ne edges as weights of CG network edges
# Max:  take the max. cij as weights of CG network edges

remodel.cna <- function(x, member = NULL, col = NULL, minus.log = TRUE,
       method = c("none", "mean", "max"), ne=3, scut=2, normalize = TRUE, 
       vmd.color = TRUE, ...) {

   require(igraph)
   method <- match.arg(method)
 
   # assume a network ensemble 
   if(inherits(x, "cna")) x = list(x)
   
   # return network names
   names <- names(x)

   # membership
   if(is.null(member)) {
      m.all = lapply(x, function(net) net$communities$membership) 
   } else {
      m.all = member
   }

   if(!is.list(m.all)) m.all = lapply(1:length(x), function(i) m.all)
 
   # number of communities 
   n <- unique(sapply(m.all, max))
   if(length(n) > 1) stop("Unequal community partition across networks")
 
   # node colors 
   if(is.null(col)) {
      if(!is.null(member)) col = 1:n
   } else {
      if(length(col) != n) 
         stop("Length of color vector doesn't match number of communities")
   }
   if(is.numeric(col) && vmd.color) col = vmd.colors()[col]

   # rebuild community network with updated cij,
   # node sizes, and edge colors
   if(!is.null(member) || method != "none") {

      col.cg.edge = NULL
      cg.node.size <- lapply(m.all, table)

      n2 <- pairwise(n)

      cg.cij <- lapply(1:length(x), function(yi) {
         y = x[[yi]]
         m = m.all[[yi]]
         cij <- y$cij
         if(minus.log) cij[cij>0] = exp(-cij[cij>0])
         if(method != "none") { 
            cij[diag.ind(cij, n=scut)] <- 0
            cij[lower.tri(cij)] <- t(cij)[lower.tri(cij)]
         }
         w <- apply(n2, 1, function(i) {
           ind1 <- which(m %in% i[1])
           ind2 <- which(m %in% i[2])
           switch(method, 
              none = max(cij[ind1, ind2]),
              max = max(cij[ind1, ind2]),
              mean= mean(sort(cij[ind1, ind2], decreasing=TRUE)[1:ne]) 
           )
         } )    
         key <- apply(n2, 1, function(i) {
           ind1 <- which(m %in% i[1])
           ind2 <- which(m %in% i[2])
           k <- switch(method, 
              none = which.max(cij[ind1, ind2]),
              max = which.max(cij[ind1, ind2]),
              mean= order(cij[ind1, ind2], decreasing=TRUE)[1:ne] 
           )
           k2 <- floor((k-1)/length(ind1)) + 1
           k1 <- (k-1) %% length(ind1) + 1
           list(i=ind1[k1], j=ind2[k2])
         } )
         cg.cij <- matrix(1, nrow=max(m), ncol=max(m))
         cg.cij[lower.tri(cg.cij)] <- w
         cg.cij[upper.tri(cg.cij)] <- t(cg.cij)[upper.tri(cg.cij)]
         list(cij = cg.cij, key = key)
      } )
      ii <- lapply(cg.cij, function(x) as.vector(sapply(x$key, "[[", "i")))
      jj <- lapply(cg.cij, function(x) as.vector(sapply(x$key, "[[", "j")))
      key <- lapply(1:length(ii), function(i) cbind(ii[[i]], jj[[i]]))
      cg.cij <- lapply(cg.cij, "[[", "cij")

      if(method != "none" && normalize) {
         w <- normalize.cij(cg.cij, ...)

         if(!is.null(w$color)) {
            col.cg.edge <- lapply(w$color, function(x) x[lower.tri(x)])
            tcol <- unlist(col.cg.edge)
            if(length(unique(tcol[!is.na(tcol)])) == 1)
               tcol[!is.na(tcol)] <- rep("#0000FF", length(sum(!is.na(tcol))))
            else 
               tcol[!is.na(tcol)] <- vec2color(tcol[!is.na(tcol)])
            tcol <- split(tcol, f=rep(1:length(col.cg.edge), each=length(col.cg.edge[[1]])) )
            col.cg.edge <- lapply(tcol, function(x) x[!is.na(x)]) 
         }
         cg.cij <- w$cij
      }

      x <- lapply(1:length(x), function(i) {
         y = x[[i]]
         if(is.list(cg.cij)) cij <- cg.cij[[i]]
         else cij <- cg.cij
         if(minus.log) {
            cij[cij>=1] <- 0.9999
            cij[cij>0] <- -log(cij[cij>0])
         }
         y$community.network <- graph.adjacency(cij, 
                                    mode = "undirected",
                                    weighted = TRUE,
                                    diag = FALSE)
         y$community.network <- set.vertex.attribute(y$community.network, "size", value=cg.node.size[[i]])
         if(!is.null(col.cg.edge)) 
            y$community.network <- set.edge.attribute(y$community.network, "color", value=col.cg.edge[[i]])
         y$communities$membership <- m.all[[i]]
         y$community.cij <- cij
         if(!is.null(key)) {
            inds = which(abs(y$cij[key[[i]]]) > 0)
            y$community.key.cij <- key[[i]][inds, ]
         }
         y
      } )
   }

   # update node color, community network node color, and community reindex
   if(!is.null(col)) {
      x <- lapply(1:length(x), function(yi) {
         y = x[[yi]]
         m = m.all[[yi]]
         y$network <- set.vertex.attribute(y$network, "color", value= col[m])
         y$community.network <- set.vertex.attribute(y$community.network, "color", value = col)
         if(!is.null(y$community.reindex)) {
            if(vmd.color) y$community.reindex = match(col, vmd.colors())
            else y$community.reindex = 1:n
         }
         y
      } )
   }

   names(x) <- names
   return(x)
}
