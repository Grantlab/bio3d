`plot.dccm` <-function(x, sse=NULL, colorkey=TRUE,
                     at=c(-1, -0.75, -0.5,  -0.25, 0.25, 0.5, 0.75, 1),
                     main="Residue Cross Correlation", # pad=0.022
                     helix.col = "gray20", sheet.col = "gray80",
                     inner.box=TRUE, outer.box=FALSE,
                     xlab="Residue No.", ylab="Residue No.",...) {

  require(lattice)
  require(grid)
  
  p1 <- contourplot(x, region = TRUE, labels=F, col="gray40",
                    at=at, xlab=xlab, ylab=ylab,
                    colorkey=colorkey,
                    main=main, ...)
  
  if (is.null(sse)) {
    print(p1)
  } else {
    xlim <- p1$x.limits
    ylim <- p1$y.limits
    uni <- 1/(max(xlim)-min(xlim))

    pad=0.02 ## This should be setable!
    xymax <- 1-(pad)#-0.002)
    pad <- pad/uni
    p1$x.limits[2]=p1$x.limits[2]+pad
    p1$y.limits[2]=p1$y.limits[2]+pad
    
    print(p1)
    
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

    
    ##- Box the plot
    if(!outer.box) {
      grid.rect(x=0, y=0,
                gp = gpar(fill=NA,col="white"),
                just=c("left","bottom"),
                width=1,height=1,
                vp=vpPath("plot_01.toplevel.vp","plot_01.panel.1.1.vp"))
    }
    if(inner.box) {
      grid.rect(x=0, y=0,
                gp = gpar(fill=NA,col="black"),
                just=c("left","bottom"),
                width=xymax, height=xymax,
                vp=vpPath("plot_01.toplevel.vp","plot_01.panel.1.1.vp"))
    }
    

  }
  ##return(p1)
}

