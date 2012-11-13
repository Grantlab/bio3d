`plot.bio3d` <-
function(x, y=NULL, type="h",
                       main="", sub="",
                       xlim=NULL, ylim=NULL, ylim2zero=TRUE,
                       xlab = NULL, ylab = NULL, 
                       axes=TRUE, ann=par("ann"),
                       col=par("col"),
                       sse=NULL, top=TRUE, bot=TRUE,
                       helix.col="gray20", sheet.col="gray80",
                       sse.border=FALSE,
                       ...) {
  
  xy <- xy.coords(x, y)
  if (is.null(xlim))
    xlim <- range(xy$x[is.finite(xy$x)])
  if (is.null(ylim))
    ylim <- range(xy$y[is.finite(xy$y)])
  if(ylim2zero) ylim[1]=0
  
  #opar <- par(no.readonly=TRUE)
  #on.exit(par(opar))
  
  plot.new()
  plot.window(xlim, ylim, ...)
  points(xy$x, xy$y, col=col, type=type, ...)

  ##- sse annotation
  if(!is.null(sse)) {
    if(top) {
      bo <- max(ylim) + (diff(ylim)*0.001) # 0.1% 
      to <- max(ylim) + (diff(ylim)*0.04) # 4%

      if(!is.null(sse$helix$start))
        rect(xleft=sse$helix$start, xright=sse$helix$end,
             ybottom=bo, ytop=to, col=helix.col, border=sse.border)#,...)

      if(!is.null(sse$sheet$start))
        rect(xleft=sse$sheet$start, xright=sse$sheet$end,
             ybottom=bo, ytop=to, col=sheet.col, border=sse.border)#,...)
    }
    if(bot){
      to <- min(ylim) - (diff(ylim)*0.001)
      bo <- min(ylim) - (diff(ylim)*0.04) # 4%

      if(!is.null(sse$helix$start))
        rect(xleft=sse$helix$start, xright=sse$helix$end,
             ybottom=bo, ytop=to, col=helix.col, border=sse.border)#,...)
      
      if(!is.null(sse$sheet$start))
        rect(xleft=sse$sheet$start, xright=sse$sheet$end,
             ybottom=bo, ytop=to, col=sheet.col, border=sse.border)#,...)
      
    }
  }
  ##- end
  if (axes) {
    axis(1)
    axis(2)
    box()
  }
  if (ann) {
    if(is.null(xlab))  xlab=xy$xlab
    if(is.null(ylab))  ylab=xy$ylab
    title(main=main, sub=sub, 
          xlab=xlab, ylab=ylab, ...)
  }
}

