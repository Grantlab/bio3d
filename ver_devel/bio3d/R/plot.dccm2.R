plot.dccm2 <-function(x, sse=NULL, colorkey=TRUE,
                     at=c(-1, -0.75, -0.5,  -0.25, 0.25, 0.5, 0.75, 1),
                     main="Residue Cross Correlation", # pad=0.022
                     helix.col = "gray20", sheet.col = "gray80",
                     inner.box=TRUE, outer.box=FALSE,
                     xlab="Residue No.", ylab="Residue No.",
                     margin.segments=NULL, segment.col=NULL, ...) {

  ## cij <- dccm(xyz)
  ## net <- igraph.comms(cij)
  ## plot.dccm2(cij, margin.segments=net$membership)
  
  require(lattice)
  require(grid)
  
  p1 <- contourplot(x, region = TRUE, labels=F, col="gray40",
                    at=at, xlab=xlab, ylab=ylab,
                    colorkey=colorkey, main=main, ...)

  xymin=0; xymax=1
  if (is.null(sse) && is.null(margin.segments)) {
    print(p1)
  } else {
    xlim <- p1$x.limits
    ylim <- p1$y.limits
    uni <- 1/(max(xlim)-min(xlim))
    pad=0.02 ## This should be setable!
    padref <- pad/uni
    
    if(!is.null(sse)) {
      ## Adjust Top and Right margins for 'sse'
      xymax <- 1-(pad)
      p1$x.limits[2]=p1$x.limits[2]+padref
      p1$y.limits[2]=p1$y.limits[2]+padref
    }
    if(!is.null(margin.segments)) {
      ## Adjust Bottom and Left margins for 'segments'
      xymin = pad
      p1$x.limits[1]=p1$x.limits[1]-padref
      p1$y.limits[1]=p1$y.limits[1]-padref

      ## Format margin annotation object
      grps <- table(margin.segments)
      small.grps=5  ## Exclude small grps
      grps = names( grps[grps > small.grps] )

      store <- NULL
      for(i in 1:length(grps)) {
        store <- rbind(store,
          cbind( bounds(which(margin.segments == grps[i])), "grp"=i))
      }
      margin.segments = store
      
      ## Margin segment colors
      if(is.null(segment.col)) {
        segment.col <- (margin.segments[,"grp"])
      } else {
        segment.col <- segment.col[(margin.segments[,"grp"])]
      }
    }
    print(p1)

    if(!is.null(sse)) {
      ##- SSE
      ##- TOP helix
      grid.rect(x=unit(sse$helix$start, "native"),
                y=xymax,
                gp = gpar(fill=helix.col,col=NA),
                just=c("left","bottom"),
                width=unit(sse$helix$length-1, "native"),
                height=unit(1, "npc"),
                vp=vpPath("plot_01.toplevel.vp","plot_01.panel.1.1.vp"))
      ## sheet
      grid.rect(x=unit(sse$sheet$start, "native"),
                y=xymax,
                gp = gpar(fill=sheet.col,col=NA),
                just=c("left","bottom"),
                width=unit(sse$sheet$length-1, "native"),
                height=unit(1, "npc"),
                vp=vpPath("plot_01.toplevel.vp","plot_01.panel.1.1.vp"))
    
      
      ##- RIGHT helix
      grid.rect(x=xymax, y=unit(sse$helix$start, "native"),
                gp = gpar(fill=helix.col,col=NA),
                just=c("left","bottom"),
                width=unit(1, "npc"),
                height=unit(sse$helix$length-1, "native"),
                vp=vpPath("plot_01.toplevel.vp","plot_01.panel.1.1.vp"))
      ## Sheet
      grid.rect(x=xymax, y=unit(sse$sheet$start, "native"),
                gp = gpar(fill=sheet.col,col=NA),
                just=c("left","bottom"),
                width=unit(1, "npc"),
                height=unit(sse$sheet$length-1, "native"),
                vp=vpPath("plot_01.toplevel.vp","plot_01.panel.1.1.vp"))
    }
    
    if(!is.null(margin.segments)) {
      ## Cluster annotation 
      ##- BOTTOM
      grid.rect(x=unit(margin.segments[,"start"], "native"),
                y=0,
                gp = gpar(fill=segment.col, col=NA),
                just=c("left","bottom"),
                width=unit(margin.segments[,"length"]-1, "native"),
                height=xymin,#height=unit(1, "npc"),
                vp=vpPath("plot_01.toplevel.vp","plot_01.panel.1.1.vp")) 

      ##- LEFT
      grid.rect(x=0, 
                y=unit(margin.segments[,"start"], "native"),
                gp = gpar(fill=segment.col, col=NA),
                just=c("left","bottom"),
                width=xymin, #width=unit(1, "npc"),
                height=unit(margin.segments[,"length"]-1, "native"),
                vp=vpPath("plot_01.toplevel.vp","plot_01.panel.1.1.vp"))
    }
    if(!outer.box) {
      grid.rect(x=0, y=0,
                gp = gpar(fill=NA,col="white"),
                just=c("left","bottom"),
                width=1,height=1,
                vp=vpPath("plot_01.toplevel.vp","plot_01.panel.1.1.vp"))
    }
    if(inner.box) {
      grid.rect(x=xymin, y=xymin,
                gp = gpar(fill=NA,col="black"),
                just=c("left","bottom"),
                width=xymax, height=xymax,
                vp=vpPath("plot_01.toplevel.vp","plot_01.panel.1.1.vp"))
    }
  }
}

