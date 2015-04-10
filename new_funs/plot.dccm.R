plot.dccm <- function(x, resno=NULL, sse=NULL, colorkey=TRUE, 
                     show = c("full", "upper", "lower"),
                     at=c(-1, -0.75, -0.5,  -0.25, 0.25, 0.5, 0.75, 1),
                     main="Residue Cross Correlation", 
                     helix.col = "gray20", sheet.col = "gray80",
                     inner.box=TRUE, outer.box=FALSE,
                     xlab="Residue No.", ylab="Residue No.",
                     margin.segments=NULL, segment.col=vmd.colors(), 
                     segment.min=1, pad=0.02, ...) {

  if(!requireNamespace("lattice", quietly=TRUE))
     stop("Please install lattice package from CRAN")

  require(grid) # should be removed after adding to the package
 
  colnames(x) = NULL; rownames(x)=NULL

  show <- match.arg(show)

  ##-- Customized axis tick labels
  if(!is.null(resno)) {
    if(is.pdb(resno)) {
      ## Take Calpha residue numbers from PDB input
      ca.inds <- atom.select(resno, "calpha", verbose = FALSE)
      resno <- resno$atom$resno[ca.inds$atom]
    }
    if(length(resno) != nrow(x)) {
      warning("Length of input 'resno' does not equal the length of input 'x'; Ignoring 'resno'")
      resno=NULL
    }
  }
  xy.at <- pretty(seq_along(x[1, ]))
  xy.at <- xy.at[xy.at <= ncol(x)]
  xy.at[1] <- 1
  if(is.null(resno)) {
     scales <- list(at=xy.at, labels=xy.at)
  } else {
     labs <- resno[xy.at]
     labs[is.na(labs)] <- ""
     scales <- list(at=xy.at, labels=labs)
  }

  ##-- Customize default plot based on show type
  if(show == "full") {
      sse.pos <- c("top", "right")
      segment.pos <- c("bottom", "left")
      diag.axis = FALSE
      if(colorkey) colorkey <- list(space="right", 
         axis.line=list(col="black"))
  } else {
     scales = c(scales, list(draw=FALSE))
     outer.box = FALSE
     diag.axis = TRUE
     diag.label = xlab
     xlab = ""
     ylab = ""
     if(show == 'upper') {
         x[row(x) > col(x)] <- NA
         sse.pos <- c("top", "left")
         segment.pos <- c("top", "left")
         if(colorkey) colorkey <- list(space="left", 
            axis.line=list(col="black"))
     } 
     if(show == 'lower') {
         x[row(x) < col(x)] <- NA
         sse.pos <- c("bottom", "right")
         segment.pos <- c("bottom", "right")
         if(colorkey) colorkey <- list(space="right", 
            axis.line=list(col="black"))
     }
  }
 
  ##-- Main Plot
  p1 <- lattice::contourplot(x, region=TRUE, labels=FALSE, col="gray40",
                    at=at, xlab=xlab, ylab=ylab,
                    colorkey=colorkey, main=main, scales=scales, ...)

  ##-- Axis on diagonal
  if(diag.axis)
     p1 <- update(p1, diag.scales=list(at=scales$at, labels=scales$labels, 
       diag.label=diag.label))

  ##-- re-write panel function for drawing segment (e.g. sse) annotations
  p1 <- update(p1, panel=.new.lattice.panel, init.limits=c(p1$x.limits, p1$y.limits))

  ##-- Customize outlines
  if(!outer.box) {
     p1 <- update(p1, par.settings = list(axis.line=list(col=NA)),
        axis = function(..., line.col) {
          lattice::axis.default(..., line.col = "black")
     } )
     if(inner.box) {
        p1 <- update(p1, inner.box = list(type=show, draw=TRUE))
     }
  } else {
     if(inner.box) {
        p1 <- update(p1, inner.box = list(type=show, draw=FALSE))
     }
  }

  ##-- Add SSE/segment annotations
  if(!is.null(sse) || !is.null(margin.segments)) {
     if(!is.null(sse))
        p1 <- add.sse(sse, pos=sse.pos, helix.col = helix.col,
                sheet.col = sheet.col, p=p1, segment.min = segment.min, pad = pad)
     if(!is.null(margin.segments))
        p1 <- add.segment(margin.segments, pos=segment.pos,
                segment.col = segment.col, p=p1, segment.min = segment.min, pad = pad)
  }

  p1
}
