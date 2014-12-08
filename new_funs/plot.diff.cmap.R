plot.diff.cmap <- function(x, y, col=c("red", "blue"), sse=NULL, sse.y=sse, col.common = "gray80", ...) {

  if(any(dim(x) != dim(y))) 
     stop("Contact maps to compare have different dimension")

  # Common contact
  comm = array(NA, dim(x))
  comm[x==1 & y==1] = 1
  comm[x==0 | y==0] = 0
 
  # Specific contact
  x.spec = x
  x.spec[x==1 & y==1] = 0
  y.spec = y
  y.spec[x==1 & y==1] = 0

  op = par(no.readonly = TRUE)
  dots = list(...)
  # reset x/y labels for this function
  if(is.null(dots$xlab)) dots$xlab = "Residue index for protein 1"
  if(is.null(dots$ylab)) dots$ylab = "Residue index for protein 2"

  # plot the common cmap; On this step, we add axes, sse (if available) and all labels.
  dots2 = dots
  if(!is.null(dots$cex)) cex = dots$cex / 2
  else cex = 0.5
  dots2$cex = cex 
  dots2$pch = 20
  # bottom sse and data, but no annotatios and axes
  dots2$left = FALSE
  dots2$axes = FALSE
  dots2$ann = FALSE
  do.call(plot.cmap, c(list(x=comm, col=col.common, sse=sse), dots2))
  # left sse and annotations/axes, but no data
  dots2$left = dots$left
  dots2$bot = FALSE
  dots2$axes = dots$axes
  dots2$ann = dots$ann
  dots2$type = "n"
  par(new=TRUE) 
  do.call(plot.cmap, c(list(x=comm, col=col.common, sse=sse.y), dots2))
  par(op)
  
  # don't repeat title, labels and x/y axes
  dots$ann=FALSE
  dots$axes=FALSE
  dots$sse=NULL

  # plot the first cmap
  par(new=TRUE) 
#  dots$add = TRUE 
  do.call(plot.cmap, c(list(x=x.spec, col=col[1]), dots))
  par(op)
  
  # plot the second cmap
  par(new=TRUE) 
#  dots$add = TRUE 
  do.call(plot.cmap, c(list(x=y.spec, col=col[2]), dots))
  par(op)

  return(invisible(NULL))
}
