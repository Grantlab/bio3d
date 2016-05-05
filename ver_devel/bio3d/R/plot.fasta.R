"plot.fasta" <- function(x, hc = TRUE, labels = x$id, cex.lab = 0.7,
                         xlab = "Alignment index",
                         main = "Sequence Alignment Overview",
                         mar4 = 4, ... ) {

  if(!(inherits(x, "fasta") | inherits(x, "pdbs")))
    stop("input 'x' should be a list object as obtained from 'read.fasta'")

  if(is.logical(hc)) {
    if(hc) {
      ide <- seqidentity(x)
      hc <- hclust(as.dist(1-ide))
    }
    else {
      hc <- NULL
    }
  }
  else {
    if(!is.null(hc)) {
      if(class(hc)!="hclust")
        stop("'hc' must be logical, NULL or a 'hclust' object.")
    }
  }
    

  if(!is.null(hc)) {
    layout(matrix(c(4,2,3,1), ncol=2),
           heights = c(.1, 1),
           widths = c(0.3, 1))
    par(mar=c(4, 0.1, 0.1, mar4))
  }
  else {
    layout(matrix(c(2, 1), nrow=2),
           heights = c(.1, 1))
    par(mar=c(4, 2, 0.1, mar4))
  }

  ## 1: gap, 0: non-gap
  gaps <- gap.inspect(x$ali)

  mat <- gaps$bin
  if(any(mat==1)) {
    mat[ mat == 1 ] <- -1
    mat[ mat == 0 ] <- 1
    mat[ mat == -1 ] <- 0
  }
  else {
    mat <- mat+1
  }

  ## re-order matrix
  if(!is.null(hc)) 
    mat <- mat[ hc$order, ]
  else
    mat <- mat[ seq(nrow(mat), 1), ]

  if(any(mat==0))
    col <- c("#FFFFFF", "#9F9F9F")
  else
    col <- c("#9F9F9F")

  image(t(mat), col = col, axes=FALSE)

  by <- pretty(0:ncol(x$ali), n = 6)
  by <- by[2]

  labs.x <- seq(0, ncol(x$ali), by=by)
  labs.x[1] <- 1

  at <- labs.x / labs.x[length(labs.x)]
  at[1] <- 0

  axis(1, at=at, labels=labs.x)
  mtext(xlab, 1, line=2.5, cex=1.0)
  at <- seq(0, 1, length.out=length(x$id))

  #labs <- x$id
  if(!is.null(hc)) 
    labels <- labels[hc$order]

  mtext(labels, side=4, line=2-1.25, at=at, cex=cex.lab, las=2)


  ## cluster dendrogram
  if(!is.null(hc)) {
    par(mar=c(4, 0.1, 0.1, 0.1))
    ddr <- as.dendrogram(hc)
    plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
  }

  ## conservation bar on top of alignment
  cons <- conserv(x$ali)
  ng <- rep(0, length(cons))
  ng[gaps$f.inds] <- 1

  if(!is.null(hc))
    par(mar=c(.1, 0.1, 2, mar4))
  else
    par(mar=c(.1, 2, 2, mar4))

  #if(all(c(0,1) %in% cons))
  #  col <- c("#FFFFFF", "#FF0000")
  #else
  col <- colorRampPalette(c("white", "red"))( 10 )

  image(as.matrix(cons), col = col, axes=FALSE)

  #if(!is.null(hc))
  #  at <- 0.4
  #else

  at <- NA
  mtext(main, side = 3, line = 0.5, cex = 1.25, at=at)
}
