plot.cna <- function(x, weights=NULL,
                     layout=NULL, col=NULL, full=FALSE, ...) {

  ##- Function for ploting cna networks the way we like them.
  ##   Returns the plot layout coordinates silently
  ##
  ## ToImprove:
  ##   Should have layout be user setable, with the options of: 
  ##  'layout.fruchterman.reingold', 'layout.mds' or 'layout.svd'
  ##
  ## \dots can contain:
  ##  col=vmd.colors(),
  ##  mark.col=vmd.colors(alpha=0.3),
  ##  mark.border=vmd.colors()
  ## AND
  ##  vertex.size:  Node sizes:   V(x$network)size
  ##  vertex.color: Node colors:  V(x$network)color
  ##  vertex.label: Node labels:  vertex ids
  ##  edge.width:   Edge weights: E(x$network)$weight
  ##  edge.color:   Edge colors:  E(x$network)$color
  ##    (also vertex.label.color, vertex.label.cex etc.
  ##     see ?igraph.plotting)
  ##
  ## col=V(x$clustered.network)$color
  ##    vmd.colors(18)== V(x$clustered.network)$color

  if(full) {
    y <- x$network 
  } else {
    y <- x$community.network
  }

  if(is.null(weights)){
    weights <- E(y)$weight
    
    if(is.null(x$call$minus.log)){  ## This stuff wont work!!!
      weights <- exp(-weights)
    }
    else{
      if(x$call$minus.log){
        weights <- exp(-weights)
      }
    }
     weights <- weights*10
  }  ## The above code is crap - should scale the weights to fit plot!!
  

  if(is.null(layout)) {
      coords <- layout.fruchterman.reingold(y, weights=weights)
  } else { coords=layout }
  
  plot.igraph(y, edge.width=weights, layout=coords, vertex.color=col, ...)
    
  ## Silently return plot coordinates
  #class(coords) = "cna"
  layout <- coords
}
