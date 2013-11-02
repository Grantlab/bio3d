plot.cna <- function(x, community=x$clustered.communities, weights=E(x$clustered.network)$weight*10,
                     layout=NULL, col=NULL, ...) {

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
  
  
  if(is.null(layout)) {
      coords <- layout.fruchterman.reingold(x$clustered.network, weights=weights)
  } else { coords=layout }
  
  if(is.null(community)) {
    ## No community structure so call plot.igraph()
    plot.igraph(x$clustered.network, edge.width=weights, layout=coords, vertex.color=col) ##, ...)
    ## Note: col=ERROR here, should be vertex.color= for plot.igraph

  } else {
    ## Plot with community groups marked 
    plot.communities(community, x$clustered.network, edge.width=weights, layout=coords,  col=col, ...)
    ###plot.communities(community, x$clustered.network, edge.width=weights, layout=coords,
    ###                 col=vmd.colors(), mark.col=vmd.colors(alpha=0.3), mark.border=vmd.colors() )
  }
  ## Silently return plot coordinates
  class(coords) = "cna"
  layout <- coords
}
