plot.diff.cmap <- function(x, y, col=c("red", "blue"), col.common = "gray80", ...) {

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

  # plot the common cmap 
  if(!is.null(dots$cex)) cex = dots$cex / 2
  else cex = 0.5
  dots2 = dots
  dots2$cex = cex 
  dots2$pch = 20
  do.call(plot.cmap, c(list(x=comm, col=col.common), dots2))

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
