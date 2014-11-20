plot.bio3d <- function(x, resno=NULL, type="h",
                  main="", sub="",
                  xlim=NULL, ylim=NULL, ylim2zero=TRUE,
                  xlab = "Residue", ylab = NULL, 
                  axes=TRUE, ann=par("ann"),
                  col=par("col"),
                  sse=NULL, sse.type="classic", sse.min.length=5,
                  top=TRUE, bot=TRUE,
                  helix.col="gray20", sheet.col="gray80",
                  sse.border=FALSE,
                  ...) {

  
  pdb2sse <- function(pdb) {
    ##- Function to obtain an SSE sequence vector from a PDB object
    ##   Result similar to that returned by stride(pdb)$sse and dssp(pdb)$sse
    ##   This could be incorporated into read.pdb() if found to be more generally useful
    
    if(is.null(pdb$helix) & is.null(pdb$sheet)) {
      warning("No helix and sheet defined in input 'sse' PDB object: try using dssp()")
      ##ss <- try(dssp(pdb)$sse)
      ## Probably best to get user to do this separately due to possible 'exefile' problems etc..
      return(NULL)
    }
    rn <- pdb$atom[pdb$calpha, c("resno", "chain")]
    ss <- rep(" ", nrow(rn))
    names(ss) = paste(rn$resno,rn$chain,sep="_")
    
    for(i in 1:length(pdb$helix$start)) {
      ss[ (rn$chain==pdb$helix$chain[i] &
           rn$resno >= pdb$helix$start[i] &
           rn$resno <= pdb$helix$end[i])] = "H"
    }
    for(i in 1:length(pdb$sheet$start)) {
      ss[ (rn$chain==pdb$sheet$chain[i] &
           rn$resno >= pdb$sheet$start[i] &
           rn$resno <= pdb$sheet$end[i])] = "E"
    }
    return(ss)
  }


  if(!is.null(resno)) {
    if(is.pdb(resno)) {
      ## Take Calpha residue numbers from PDB input
      ca.inds <- atom.select(resno, "calpha", verbose = FALSE)
      resno <- resno$atom$resno[ca.inds$atom]
    }
    if(length(resno) != length(x)) {
      warning("Length of input 'resno' does not equal the length of input 'x'; Ignoring 'resno'")
      resno=NULL
    }
  }
  
  xy <- xy.coords(x)
  if (is.null(xlim))
    xlim <- range(xy$x[is.finite(xy$x)])
  if (is.null(ylim))
    ylim <- range(xy$y[is.finite(xy$y)])
  if(ylim2zero) ylim[1]=0
  
  plot.new()
  plot.window(xlim, ylim, ...)
  points(xy$x, xy$y, col=col, type=type, ...)
  
  if(!is.null(sse)) {	

    ## Obtain SSE vector from PDB input 'sse'
    if(is.pdb(sse)) 
      sse$sse <- pdb2sse(sse)
 
    h <- bounds( which(sse$sse == "H") )
    e <- bounds( which(sse$sse == "E") )
    
    ## Remove short h and e elements that can crowd plots
    if(length(h) > 0) {
      h <- h[h[,"length"] >= sse.min.length,]
    } else { h <- NULL }
    if(length(e) > 0) {
      e <- e[e[,"length"] >= sse.min.length,]
    } else { e <- NULL }

    if(sse.type != "classic")
      warning("Only sse.type='classic' is currently available, 'fancy' coming soon")
    
    if(top) {
      ## Determine bottom and top of margin region 
      bo <- max(ylim) + (diff(ylim)*0.001) # 0.1% 
      to <- max(ylim) + (diff(ylim)*0.04) # 4%
      
      if(length(h) > 0)
        rect(xleft=h[,"start"], xright=h[,"end"],
             ybottom=bo, ytop=to, col=helix.col, border=sse.border)
      
      if(length(e) > 0)
        rect(xleft=e[,"start"], xright=e[,"end"],
             ybottom=bo, ytop=to, col=sheet.col, border=sse.border)
    }
    if(bot){
      	to <- min(ylim) - (diff(ylim)*0.001)
      	bo <- min(ylim) - (diff(ylim)*0.04)
        
      	if(length(h) > 0)
          rect(xleft=h[,"start"], xright=h[,"end"],
               ybottom=bo, ytop=to, col=helix.col, border=sse.border)
        
      	if(length(e) > 0)
          rect(xleft=e[,"start"], xright=e[,"end"],
               ybottom=bo, ytop=to, col=sheet.col, border=sse.border)
      }
  }

  if(axes) {
    axis(2)
    box()
    at <- axTicks(1); at[1] = 1
    if(is.null(resno)) {
      axis(1, at)
    } else {
      labels <- resno[at]
      axis(1, at=at, labels=labels)
    }
  }
  if(ann) {
    if(is.null(xlab))  xlab=xy$xlab
    if(is.null(ylab))  ylab=xy$ylab
    title(main=main, sub=sub, 
          xlab=xlab, ylab=ylab, ...)
  }
}


