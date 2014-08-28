"plot.hclust2" <- function(hc, k=NULL, h=NULL, labels=NULL, columnColors=TRUE, ...) {

  if(!inherits(hc, "hclust"))
    stop("hc must be of type 'hclust'")
  if(is.null(k) & is.null(h))
    stop("provide either k or h to function 'cutree'")
  
  mtext.names <- names(formals( mtext ))
  plot.names <- c(names(formals( graphics:::plotHclust )),
                  names(formals( plot.default )))
  
  dots <- list(...)
  mtext.args <- dots[names(dots) %in% mtext.names]
  plot.args <- dots[names(dots) %in% plot.names]
  par.args <- dots[!(names(dots) %in% unique(c(names(mtext.args), names(plot.args))))]

  mtext.args <- c(mtext.args, par.args)
  plot.args <- c(plot.args, par.args)

  if(!any(names(plot.args)=="hang"))
    plot.args$hang <- -1
  if(!any(names(mtext.args)=="line")) {
    if( columnColors )
      mtext.args$line <- 0
    else
      mtext.args$line <- -1
  }
  if(!any(names(mtext.args)=="side"))
    mtext.args$side <- 1
  if(!any(names(mtext.args)=="las"))
    mtext.args$las <- 2
  if(any(names(mtext.args)=="col"))
    mtext.args$col <- NULL
  if(any(names(mtext.args)=="at"))
    mtext.args$at <- NULL
  
  ## print(mtext.args)
  ## print(par.args)
  ## print(plot.args)

  plot.labels <- TRUE
  if(is.logical(labels)) {
    if(labels)
      plot.labels <- TRUE
    else
      plot.labels <- FALSE
    labels <- NULL
  }
  
  if(is.null(labels)) {
    labels <- hc$labels
    if(is.null(labels))
      labels <- seq(1, length(hc$order))
  }
  else {
    if( length(hc$order) != length(labels) )
      stop("labels must be of same length as hc")
  }
  
  grps        <- cutree(hc, k=k, h=h)
  hcd         <- as.dendrogram(hc)
  labelColors <- seq(1, length(unique(grps)))
  cols <- labelColors[grps][hc$order]
  
  do.call('plot',  c(list(x=hc, labels=FALSE), plot.args))

  if(plot.labels) {
    do.call('mtext', c(list(text=labels[ hc$order ], 
                            at=1:length(grps),
                            col=cols),
                     mtext.args))
  }

  if( columnColors ) {
    for(i in 1:length(unique(cols))) {
      inds <- which(cols==i)
      rect(inds[1]-.5, -.2, inds[length(inds)]+.5, -.05, col=i)
    }
  }
  
  ##do.call('mtext', c(list(text=labels, side = 1, line=0,
  ##                        col=labelColors[clusMember][hc$order],
  ##                        at=1:length(clusMember), las=2),
  ##                   mtext.args))

}
