print.cna <- function(x, ...) {

  y <- summary.cna(x, ...)

  l1 <- paste( "\nRAW NODES#:", x$communities$vcount,
              "   COMMUNITIES#:", max(x$communities$membership),
              "   EDGES#:", ecount(x$network))


  l2 <- paste("   EDGES#:", ecount(x$community.network))

  
  i <- paste( attributes(x)$names, collapse=", ")
  cat(strwrap(paste(" + attr:",i,"\n"),width=60, exdent=8), sep="\n")
  cat(l1,"\n")
  print.igraph(x$network)
  cat(l2,"\n")
  print.igraph(x$community.network)

}
