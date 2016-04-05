col2hex <- function (cname) 
{
    colMat <- col2rgb(cname)
    rgb(red = colMat[1, ]/255, green = colMat[2, ]/255, blue = colMat[3, 
        ]/255)
}
vec2color <- function(vec, pal=c("blue", "green", "red"), n) {
  ##-- Define a color scale from a numeric vector
  require(classInt)
  return( findColours(classIntervals(vec, n=n, style="equal"), pal) )
}

get.edgecolor <- function(cij, method, colmap, cutoff, n = 30, signif = NULL, p.cutoff=0.05, ...) {

   if(length(cij) == 1 || method == 'none') return(NULL)

   require(classInt)

   pcij <- do.call(rbind, lapply(cij, function(x) x[lower.tri(x)]))

   if(!is.null(signif)) {
      signif <- signif[lower.tri(signif)]
      signif[is.nan(signif)] <- 1
   }
   color.code <- switch(method,  
      variance = {
         vars = apply(pcij, 2, sd)
         rep(vars, each=nrow(pcij))
      }, 
      feature = {
         as.vector( apply(pcij, 2, function(x) {
            if(all(x==0) || ((max(x) - min(x))/max(x)) < cutoff)
               return(rep(0, length(x)))
            class <- suppressWarnings(classIntervals(x, n=2, style='equal'))
            flag <- as.integer(x >= class$brks[2])
            flag <- flag * (1:length(flag))
         }) )
      },
      'significance' = {
         as.vector( sapply(1:ncol(pcij), function(i) {
            x <- pcij[, i]
            if(all(x==0) || signif[i] > p.cutoff)
               return(rep(0, length(x)))
            class <- suppressWarnings(classIntervals(x, n=2, style='equal'))
            flag <- as.integer(x >= class$brks[2])
            flag <- flag * (1:length(flag))
         }) )
      } )
  
   colors <- rep(NA, length(color.code)) 
   color.code <- color.code[as.vector(pcij > 0)] # exclude pairs that have no edge 
   if(length(unique(color.code)) == 1)
      colors[as.vector(pcij > 0)] <- rep(colmap[1], length(color.code))
   else 
      colors[as.vector(pcij>0)] <- switch(method, 
         variance = suppressWarnings( vec2color(color.code, pal = colmap, n = n) ),
         feature = suppressWarnings( vec2color(color.code, pal = colmap[sort(unique(color.code)+1)], 
                       n = length(unique(color.code))) ),
         significance = suppressWarnings( vec2color(color.code, pal = colmap[sort(unique(color.code)+1)],
                       n = length(unique(color.code))) )
      )
   colors <- split(colors, f = rep(1:nrow(pcij), ncol(pcij)))
   colors <- lapply(colors, function(x) x[!is.na(x)])
   return(colors)
}

normalize.cij <- function(cij, ...) {
   pcij <- do.call(rbind, lapply(cij, function(x) x[lower.tri(x)]))
   pcij[pcij <= 0] <- NA
   for(i in 1:nrow(pcij)) {
      vmin <- min(pcij[i, ], na.rm = TRUE)
      vmax <- max(pcij[i, ], na.rm = TRUE)
      if(vmax > vmin)
         pcij[i, ] <- (pcij[i, ] - vmin) / (vmax - vmin) * 0.9 + 0.05
      else
         pcij[i, !is.na(pcij[i, ])] <- 0.95
   }
   pcij[is.na(pcij)] <- 0

   ncij <- apply(pcij, 1, function(x) {
      xx = cij[[1]]
      xx[lower.tri(xx)] <- x
      xx[upper.tri(xx)] <- t(xx)[upper.tri(xx)]
      return(list(xx))
   } )
   do.call("c", ncij)
}

# Methods to calculate community cijs:
#    - Mean, use the mean of top 'ne' inter-community cijs as the cij of the CG community network
#    - Max,  use the maximum inter-community cij as the cij of the community network
#    - Sum,  use the sum of all inter-community cijs as the cij of the community network
remodel.cna <- function(x,  member = NULL, col = NULL, minus.log = TRUE,
       method = c('none', 'sum', 'mean', 'max'), ne=3, scut=4, normalize = TRUE, 
       vmd.color = TRUE, col.edge=c('none', 'variance', 'feature', 'significance'),
       colmap.edge = NULL,  coledge.cutoff = 0.5, signif=NULL, cijs = NULL, one.net=FALSE, ncore=NULL, ...) {

   require(igraph)
   require(parallel)
   method <- match.arg(method)
   col.edge <- match.arg(col.edge)

   ncore = setup.ncore(ncore)

   # assume input 'x' is a network ensemble 
   if(inherits(x, "cna")) x = list(x)
   
   # will return with network names
   net.names <- names(x)

   if(col.edge == 'significance') { 
      if(length(x) == 1)
         stop('Edges cannot be colored by "significance" provided only one network')
      if(is.null(signif)) {
         if(is.null(cijs)) 
            stop('Provide either pre-calculated significance matrices or the raw residue cross-correlation matrices.')
         
         filters <- lapply(x, function(x) x$cij > 0)
         flag <- rep(1:length(x), times=sapply(cijs, function(x) dim(x)[3L]))
         cijs <- lapply(cijs, function(x) lapply(1:dim(x)[3L], function(i) x[,,i]))
         cijs <- do.call('c', cijs)
#         flag <- rep(1:length(x), each = length(cijs)/length(x))
         cg.cij <- mclapply(1:length(cijs), function(i) {
            net <- cna(cijs[[i]]*filters[[flag[i]]], cutoff.cij=0)
            nnet <- remodel.cna(net, member=member, method='sum', scut=4, normalize=FALSE, col.edge='none')[[1]]
            nnet$community.cij[lower.tri(nnet$community.cij)]
         }, mc.cores = ncore )
          
         signif <- unlist( mclapply(1:length(cg.cij[[1]]), function(i) {
            dat <- sapply(cg.cij, '[', i)
            n2 <- pairwise(length(x))
            p <- sapply(1:nrow(n2), function(j) {
               x <- dat[flag == n2[j, 1]]
               y <- dat[flag == n2[j, 2]]
               if(any(x > 0) || any(y>0))
                  t.test(x, y)$p.value
         #         wilcox.test(x, y)$p.value
         #      else if(any(x > 0)) 
         #         t.test(x, alternative='greater')$p.value
         #      else if(any(y > 0))
         #         t.test(y, alternative='greater')$p.value
               else 1
            })
            min(p)
         }, mc.cores = ncore) )
         a <- matrix(1, nrow=length(unique(member)), ncol=length(unique(member)))
         a[lower.tri(a)] <- signif
         signif=a
      }
   }
 
   # if not provided, the membership is extracted from the network
   if(is.null(member)) {
      member = lapply(x, function(y) y$communities$membership)
   } else {
      member = rep(list(member), length(x))
      if(method == 'none') {
         warning('method is "none" but member is not NULL. Set method to "sum"')
         method = 'sum'
      }
   }
 
   if(method != 'none') {
      # calculate cijs for CG networks with defined membership
      cg.cij <- lapply(1:length(x), function(i) {
         cij <- x[[i]]$cij
         member <- member[[i]]
         n.mem <- length(unique(member))
         n2 <- pairwise(n.mem)
 
         if(minus.log) cij[cij>0] = exp(-cij[cij>0])
         cij[diag.ind(cij, n=scut)] <- 0
         cij[lower.tri(cij)] <- t(cij)[lower.tri(cij)]
   
         w <- apply(n2, 1, function(i) {
           ind1 <- which(member %in% i[1])
           ind2 <- which(member %in% i[2])
           switch(method, 
              max = max(cij[ind1, ind2]),
              mean= mean(sort(cij[ind1, ind2], decreasing=TRUE)[1:ne]),  
              sum = sum(cij[ind1, ind2]) 
           )
         } )
         cg.cij <- matrix(1, nrow=n.mem, ncol=n.mem)
         cg.cij[lower.tri(cg.cij)] <- w
         cg.cij[upper.tri(cg.cij)] <- t(cg.cij)[upper.tri(cg.cij)]
         cg.cij
      } )
   } else {
      cg.cij <- lapply(x, function(y) {
         cij <- y$community.cij
         if(minus.log) cij[cij>0] = exp(-cij[cij>0])
         cij
      })
   }

   if(method != 'none' || 
      any(sapply(x, function(y) is.null(y$community.key.cij))) ) {
      ###############################################################
      # Get the indices of cijs that contribute to the community cijs.
      # This will be stored in the returned networks and will be used
      # with another application related to plot 3D networks with VMD.
      # For 2D plot, this variable is not used.
      key <- lapply(1:length(x), function(i) {
         cij <- x[[i]]$cij
         member <- member[[i]]
         n.mem <- length(unique(member))
         n2 <- pairwise(n.mem)

         if(minus.log) cij[cij>0] = exp(-cij[cij>0])
         cij[diag.ind(cij, n=scut)] <- 0
         cij[lower.tri(cij)] <- t(cij)[lower.tri(cij)]
   
         val <- apply(n2, 1, function(i) {
           ind1 <- which(member %in% i[1])
           ind2 <- which(member %in% i[2])
           k <- switch(method, 
              max = which.max(cij[ind1, ind2]),
              mean= order(cij[ind1, ind2], decreasing=TRUE)[1:ne],
              sum = order(cij[ind1, ind2], decreasing=TRUE)
           )
           k2 <- floor((k-1)/length(ind1)) + 1
           k1 <- (k-1) %% length(ind1) + 1
           list(i=ind1[k1], j=ind2[k2])
         })
         cbind(unlist(lapply(val, '[[', 'i')), unlist(lapply(val, '[[', 'j')))
      })
      #################################################################
   } else {
      key <- lapply(x, '[[', 'community.key.cij')
   }

   if(length(x) == 1) {
      edge.color = NULL
   } else {
      if(col.edge == 'none') {
         edge.color = NULL
      } else {
         check <- TRUE
         mlist <- lapply(member, function(x) sort(unique(x)))
         if(length(unique(sapply(mlist, length))) > 1) {
           check = FALSE
         } else {
           mmat <- do.call(rbind, mlist)
           if(any(apply(mmat, 2, function(x) length(unique(x)) > 1))) 
             check = FALSE
         }
         
#         if(length(unique(sapply(member, length))) != 1) check = FALSE
#         else 
#            check <- all(apply(do.call(rbind, member), 2, function(x) length(unique(x))==1))
         if(!check) edge.color=NULL
         else {
            # set default colormap according to the color method for edges
            if(is.null(colmap.edge)) {
               colmap.edge <- switch(col.edge, 
                  variance = c('blue', 'red'),
                  feature = {
                     tcol = c('gray', 'red', 'darkgreen', 'blue') 
                     tcol = union(tcol, colors()[-1])
                     tcol[1:(length(x)+1)]
                  },
                  significance = {
                     tcol = c('gray', 'red', 'darkgreen', 'blue') 
                     tcol = union(tcol, colors()[-1])
                     tcol[1:(length(x)+1)]
               })
            }
            if((col.edge %in% c('feature', 'significance')) && length(colmap.edge) != (length(x)+1))
               stop('Number of colors does not match input number of networks')

            # Calculate the variance of CG cijs across networks.
            # The values will be used to color the edges of CG networks
            edge.color <- get.edgecolor(cij = cg.cij, method = col.edge,
                   colmap = colmap.edge, cutoff = coledge.cutoff, signif=signif, ...)
         }
      }
   }
   if(is.null(edge.color)) 
      edge.color <- lapply(x, function(y) 
          get.edge.attribute(y$community.network, name='color') )

   # Normalize CG cijs for each network
   if(normalize) cg.cij <- normalize.cij(cg.cij, ...)

   # calculate node sizes and colors
   cg.node.size <- lapply(member, table)
   n.mem <- sapply(member, function(y) length(unique(y)))
   if(is.null(col)) {
      check <- all(sapply(1:length(n.mem), 
                  function(i) n.mem[i] == length(V(x[[i]]$community.network))))
      if(check) 
         col <- lapply(x, function(y) get.vertex.attribute(y$community.network, 'color'))
      else 
         col = lapply(member, function(y) {
                    col <- 1:length(unique(y))
                    if(vmd.color) vmd.colors()[col]
                    else col
               } )
   } else {
      if(!all(sapply(n.mem, '==', length(col))))
         stop("Length of color vector doesn't match number of communities")
      if(is.numeric(col) && vmd.color) col = vmd.colors()[col]
      col = rep(list(col), length(member))
   }
   # update network components
   x <- lapply(1:length(x), function(i) {
      y = x[[i]]
      cij <- cg.cij[[i]]
      if(minus.log && (method != 'sum' || normalize)) {
         cij[cij>=1] <- 0.9999
         cij[cij>0] <- -log(cij[cij>0])
      }
      y$community.network <- graph.adjacency(cij, 
                                 mode = "undirected",
                                 weighted = TRUE,
                                 diag = FALSE)
      y$community.network <- set.vertex.attribute(y$community.network, "size", value=cg.node.size[[i]])

      y$community.network <- set.edge.attribute(y$community.network, "color", value=edge.color[[i]])
      y$communities$membership <- member[[i]]
      y$community.cij <- cij
      inds = which(abs(y$cij[key[[i]]]) > 0 & abs(apply(key[[i]], 1, diff)) >= scut)
      y$community.key.cij <- cbind(key[[i]][inds, ,drop=FALSE], member[[i]][key[[i]][inds,1]],
                                                     member[[i]][key[[i]][inds,2]])

      y$network <- set.vertex.attribute(y$network, "color", value= col[[i]][member[[i]]])
      y$community.network <- set.vertex.attribute(y$community.network, "color", value = col[[i]])
      if(!is.null(y$community.reindex)) {
         if(vmd.color) y$community.reindex = match(col[[i]], vmd.colors())
         else y$community.reindex = 1:length(unique(member[[i]]))
      }
      y
   } )

   names(x) <- net.names

   if(one.net && length(x) == 2) {
     y <- x[[1]]
     size <- get.vertex.attribute(y$community.network, "size")
     col <-  get.vertex.attribute(y$community.network, "color")
     cij <- y$community.cij
     cij[y$community.cij == 0] <- x[[2]]$community.cij[y$community.cij == 0]
     y$community.network <- graph.adjacency(cij,
                           mode = "undirected",
                           weighted = TRUE,
                           diag = FALSE)

     y$community.network <- set.vertex.attribute(y$community.network, "size", value=size)
     y$community.network <- set.vertex.attribute(y$community.network, "color", value = col)
     y$community.cij <- cij

     common.color <- col2hex(colmap.edge[1])
     cij1 <- x[[1]]$community.cij
     cij2 <- x[[2]]$community.cij
     cij1 <- cij1[lower.tri(cij1)]
     cij2 <- cij2[lower.tri(cij2)]
     ecol1 <- get.edge.attribute(x[[1]]$community.network, "color")
     ecol2 <- get.edge.attribute(x[[2]]$community.network, "color")
     ecol <- rep(NA, length(cij1))
     ecol[cij1 > 0] <- ecol1
     ecol[cij2 > 0][ecol2 != common.color] <- ecol2[ecol2 != common.color]
     ecol <- ecol[!is.na(ecol)]
     y$community.network <- set.edge.attribute(y$community.network, "color", value=ecol)
     
     elty <- rep(NA, length(cij1))
     elty[cij1 > 0] <- 1
     elty[cij2 > 0 & cij1 == 0] <- 3
     elty <- elty[!is.na(elty)]
     y$community.network <- set.edge.attribute(y$community.network, "lty", value=elty)
  
     w1 <- rep(0, length(cij1))
     w2 <- rep(0, length(cij1))
     w1[cij1 > 0] <- get.edge.attribute(x[[1]]$community.network, "weight")
     w2[cij2 > 0] <- get.edge.attribute(x[[2]]$community.network, "weight")
     deltaW <- abs(w1 - w2)
     deltaW <- deltaW[cij1 > 0 | cij2 > 0]
     deltaW[ecol == common.color] <- 0 
     y$delta.community.weight <- deltaW
     x <- y
   } 

   return(x)
}
