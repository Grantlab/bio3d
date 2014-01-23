"plot.fasta" <- function(x, plot.labels=TRUE, plot.conserv=TRUE,
                         conserv.col=NULL, conserv.height.scale=1,
                         conserv.intervals=c(1, .75, .5, .25),
                         row.height=100, row.spacer=20,
                         cex.text=1, add=FALSE, ... ) {

  if(class(x)!="fasta")
    stop("input 'x' should be a list object as obtained from 'read.fasta'")

  ## calculate sequence entropy
  ent <- entropy(x)
  hn <- ent$H.norm

  ## check for gaps
  gaps <- gap.inspect(x$ali)
  dims <- dim(x$ali)
  
  
  ## set x- and y-lim
  if(plot.labels)
    xlim <- c(0, dims[2]+100)
  else
    xlim <- c(0, dims[2])

  if(plot.conserv) {
    conserv.height <- row.height * conserv.height.scale
    ylim <- c(0, dims[1]*row.height + conserv.height)
  }
  else {
    conserv.height <- NULL
    ylim <- c(0, dims[1]*row.height)
  }

  if(is.null(conserv.col))
    conserv.col <- heat.colors(length(conserv.intervals)-1)

  ## start the plot
  if(!add) {
    plot.new()
    plot.window(xlim, ylim)
  }
  
  ## plot the sequence 'boxes'
  for ( i in 1:nrow(x$ali) ) {
    nongap.inds <- which(gaps$bin[i,]==0)
    bs <- bounds(nongap.inds)
    
    for ( j in 1:nrow(bs) ) {
      xstart <- bs[j,"start"]
      xend <- bs[j,"end"]
      ystart <- (i-1)*row.height
      rect(xstart, ystart+1, xend, ystart+(row.height-row.spacer), col="grey50", border=NA)
    }

    ## label on the side of the sequence
    if(plot.labels)
      text(dims[2] + 10, ystart+(row.height*.5), labels=x$id[i], pos=4,
           cex=cex.text)
  }


  ## plot conservation bar
  if(plot.conserv) {
    i <- i+1
    
    for ( j in 1:(length(conserv.intervals)-1) ) {
      cons.inds <- intersect(
        intersect(which(hn <= conserv.intervals[j]),
                  which(hn > conserv.intervals[j+1])),      
        gaps$f.inds)
      
      bs <- bounds(cons.inds)
      
      for ( k in 1:nrow(bs) ) {
        xstart <- bs[k,"start"]
        xend <- bs[k,"end"]+1
        ystart <- (i-1)*row.height
        rect(xstart, ystart, xend, ystart+(conserv.height*conserv.intervals[j]),
             col=conserv.col[j], border=NA)
      }
      
    }

    if(plot.labels)
      text(dims[2] + 10, ystart+(row.height*.5), labels="Conservation", pos=4,
           cex=cex.text)
  }
    
}
