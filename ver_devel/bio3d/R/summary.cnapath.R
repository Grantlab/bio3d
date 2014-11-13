summary.cnapath <- function(pa, ..., pdb = NULL, label = NULL, col = NULL, 
   plot = FALSE, concise = FALSE, cutoff = 0.1, normalize = TRUE) {

   pa <- list(pa, ...)
   if(!all(sapply(pa, inherits, "cnapath")))
      stop("Input pa is not a 'cnapath' object")
   
   if(is.null(label)) label = 1:length(pa)
   if(is.null(col)) col = 1:length(pa)

   out <- list()
   
   # read node numbers on paths
   y <- lapply(pa, function(x) unlist(x$path))
   
   # store node degeneracy 
   node.deg <- lapply(y, table)
   if(normalize) {
      node.deg <- lapply(node.deg, function(x) x/max(x))
   }

   # find on-path node by the cutoff
   yy <- lapply(node.deg, function(x) x[x >= cutoff])
   onpath.node <- unique(names(unlist(yy)))
   i <- as.numeric(onpath.node)
   onpath.node <- onpath.node[order(i)]

   # generate the node degeneracy table
   o <- lapply(node.deg, function(x) {
      x <- x[match(onpath.node, names(x))]
      x[is.na(x)] <- 0
      names(x) <- onpath.node 
      x 
   } )

   # replace node id with pdb resid and resno
   if(!is.null(pdb)) {
      ca.inds <- atom.select(pdb, elety="CA", verbose = FALSE)
      resno <- pdb$atom[ca.inds$atom, "resno"]
      resid <- pdb$atom[ca.inds$atom, "resid"]
      lig.inds <- atom.select(pdb, "ligand", verbose = FALSE)
      islig <- resno %in% pdb$atom[lig.inds$atom, "resno"]
      resid[!islig] <- aa321(resid[!islig])
      o <- lapply(o, function(x) {
            n <- paste(resid[as.numeric(names(x))], resno[as.numeric(names(x))], sep="")
            names(x) <- n
            x 
      } )
   }

   names(o) <- label
   out$network <- label
   out$num.paths <- sapply(pa, function(x) length(x$path))
   out$hist <- lapply(pa, function(x) table(cut(x$dist, breaks=5, include.lowest = TRUE)))
   if(length(out$hist)==1) out$hist = out$hist[[1]]
   out$degeneracy <- do.call(rbind, o)
   if(normalize) out$degeneracy <- round(out$degeneracy, digits=2)
   
   if(plot) {
      if(!requireNamespace("ggplot2", quietly = TRUE))
         stop("Please install the ggplot2 package from CRAN")
  
      gcol = col; names(gcol) = label
      
      ##- for path length distribution
      y1 <- lapply(pa, "[[", "dist")
      factor1 <- rep(label, sapply(y1, length))
      bw1 = diff(range(unlist(y1))) * 0.02
      wd1 = bw1 / 1.5
      df1 <- data.frame(State=as.factor(factor1), Length=unlist(y1))
      df1 <- cbind(df1, bw = bw1)

      ##- for node degeneracy 
      y2 <- lapply(pa, function(x) unlist(x$path))
      if(!is.null(pdb)) y2 <- lapply(y2, function(x) resno[x])
      if(concise) { 
         # re-number node to get more concise plot
         ii <- sort(unique(unlist(y2)))
         y2 <- lapply(y2, match, ii)
      }
      factor2 <- rep(label, sapply(y2, length))
      bw2 = 0.9  # must be less than 1
      wd2 = bw2 / 1.5
      df2 <- data.frame(State=as.factor(factor2), Node=unlist(y2))
      df2 <- cbind(df2, bw = bw2)
   
      p1 <- ggplot2::ggplot(df1, ggplot2::aes(x=Length, group=State, color=NULL, fill=State, bw=bw)) + 
        ggplot2::geom_histogram(ggplot2::aes(y=(..density..)*bw), alpha=0.8, 
        position = ggplot2::position_dodge(width=wd1), binwidth=bw1) + 
        ggplot2::scale_fill_manual(values=gcol) + 
        ggplot2::ylab("Probability") + 
        ggplot2::ggtitle("Path length distribution") + ggplot2::theme_bw()

      p2 <- ggplot2::ggplot(df2, ggplot2::aes(x=Node, group=State, color=NULL, fill=State, bw=bw)) + 
         ggplot2::geom_histogram(ggplot2::aes(y=..count..), alpha=0.8, 
            position=ggplot2::position_dodge(width=wd2), binwidth=bw2) + 
               ggplot2::scale_fill_manual(values=gcol)+
               ggplot2::ylab("Degeneracy") + 
               ggplot2::ggtitle("Node degeneracy") + ggplot2::theme_bw()
  
      pushViewport(viewport(layout = grid.layout(1, 2)) )
      print(p1, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
      print(p2, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
   } 

   return(out)
}

print.cnapath <- function(pa, ...) {
   dots = list(...)

   if(is.list(pa) && all(sapply(pa, inherits, "cnapath"))) {
      if(!"label" %in% names(dots) || is.null(dots$label)) dots$label = names(pa)
      names(pa) <- NULL
      args = c(pa, dots)
      o <- do.call(summary, args)
   } else {
      o <- summary(pa, ...)
   }

   if("plot" %in% names(dots)) plot = dots$plot
   else plot = FALSE

   if(!plot) {
      if("normalize" %in% names(dots)) normalize = dots$normalize
      else normalize = TRUE
   
      if(length(o$network) > 1) {   
         cat("Number of networks: ", length(o$network), "(", 
             paste(o$network, collapse=", "), ")\n")
      }
   
      cat("Number of paths in network(s):\n")
      if(length(o$network) > 1) {
          cat(paste("   ", o$network, ": ", o$num.paths, sep="", collapse="\n"), sep="\n")
          cat("\n")
      } else {
          cat("    ", o$num.paths, "\n\n")
      }
   
      cat("Path length distribution: \n")
      if(length(o$network) > 1) {   
         for(i in 1:length(o$network)) {
             cat("   --- ", o$network[i], " ---")
             print(o$hist[[i]])
             cat("\n")
         }
      } else {
         print(o$hist)
         cat("\n")
      }
   
      cat("Node degeneracy table: \n\n")
      if(length(o$network) == 1) rownames(o$degeneracy) = ""
      if(normalize)
         print(format(o$degeneracy, nsmall=2), quote=FALSE)
      else 
         print(o$degeneracy)
   }
}
