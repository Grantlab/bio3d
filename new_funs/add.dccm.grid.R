
add.dccm.grid <- function(x, fill.col=NA, helix.col="purple", sheet.col="yellow",
                          line.col="black", lty=1, lwd=1, segment.min=1,
                          alpha=0.3, side=c(1,2)) {
  
  ## Add a grid or colored boxes to a plot.dccm() plot
  ##   sse <- dssp(pdb)
  ##   plot.dccm(cij, sse=see)
  ##   add.dccm.grid(sse)
  ##   add.dccm.grid(list("start"=c(50, 100, 150), "length"=c(10,25,50)))
  ##   add.dccm.grid(list("start"=c(50, 100, 150), "length"=c(10,25,50)),
  ##                 fill.col=c("blue","red","green"))
  ##   plot.dccm2(cij, margin.segments=net$membership, segment.min=15)
  ##   add.dccm.grid( net$membership, segment.min=15)
  ##
  ##   add.dccm.grid( net$membership, segment.min=25)
  ##   add.dccm.grid( net$membership, segment.min=25, fill.col="gray")

  ##-- Function to draw box on plot
  draw.box <- function(start, length, xymin=0, xymax=1,
                       fill.col="gray", alpha=0.3, line.col="black", lty=1, lwd=1, 
                       side=c(1,2) ) {
    
    ##-- Draw Annotation Blocks On DCCM Plots
    ##    draw.box(150,20) 
    ##    draw.box(50,100, side=1, fill.col=NA, lty=2)

    ## Grid graphics paramaters
    gp <- gpar(fill=fill.col, col=line.col,
               lty=lty, lwd=lwd, alpha=alpha)

    vp <- vpPath("plot_01.toplevel.vp",
                 "plot_01.panel.1.1.vp")

    ##- Side 1: From Bottom Margin
    if( (side==1) || (side=="both") || all(side==c(1,2)) ) {
      grid.rect(x=unit(start-0.5, "native"),
                y=xymin,
                width=unit(length-0.5, "native"),
                height=xymax,
                gp=gp, just=c("left","bottom"), vp=vp) 
    }
    
    ##- Side 2: From Left Margin
    if( (side==2) || (side=="both") || all(side==c(1,2)) ) {
      grid.rect(x=xymin, 
                y=unit(start-0.5, "native"),
                width=xymax,
                height=unit(length-0.5, "native"),
                gp=gp, just=c("left","bottom"), vp=vp)      
    }
  }


  ##NOTE: For all 'x' objects that are not vectors we will exclude
  ##      segments that are under 'segment.min' in length
  segment.min.exclusion = TRUE 
  


  ##-- Parse 'x' Input --##
  
  ##- For vector input objects - e.g. $membership from cutree()
  if( (is.vector(x)) && (!is.list(x)) ) {
    grps <- table(x)

    ## Exclude small grps less than 'segment.min'
    ##  But do not do more filtering below 
    grps = names( grps[grps > segment.min] )
    segment.min.exclusion=FALSE ##<--- good idea but plots are too crowded!!

    store.grps <- NULL; 
    for(i in 1:length(grps)) {
      store.grps <- rbind(store.grps,
          cbind( bounds(which(x == grps[i])),
                "grp"=as.numeric(grps[i])) )
    }
    ## convert to matrix for use below
    x=store.grps

    ## Dont do any more filtering
    
  }

  ##- For SSE objects
  if(class(x) == "dssp") {
    start <- c(x$helix$start, x$sheet$start)
    length <- c(x$helix$length, x$sheet$length)
    
    ## If no 'fill.col' is provided use helix and sheet specific fill.col
    if(is.na(fill.col)) {
      fill.col <- c(rep(helix.col, length(x$helix$start)),
                    rep(sheet.col, length(x$sheet$start)) )
    }
  }  

  ##- For other list objects
  if(class(x) == "list") {
    start <- x$start
    length <- x$length
  }

  ##- For matrix objects - e.g. from bounds()
  if(class(x) =="matrix") {
    start  <- x[,"start"]
    length <- x[,"length"]
  }


  ##-- Filter out short segments based on input 'segment.min'
  if(segment.min.exclusion) {
    inds <- !(length < segment.min)
    start <- start[inds]
    length <- length[inds]
    if(length(fill.col > 1)) { fill.col <- fill.col[inds] }
  }
  
  ## Draw
  draw.box(start, length, fill.col=fill.col, alpha=alpha)
}

