print.cna <- function(x, ...) {

  ## ToDo: !! Should also print details of x$clustered.network !!
  
  y <- summary.cna(x, ...)

  l1 <- paste( "\nRAW NODES#:", x$raw.communities$vcount,
              "   COMMUNITIES#:", max(x$raw.communities$membership),
              "   EDGES#:", ecount(x$raw.network))

  l2 <- paste( "\nCLUSTERED NODES#:", x$clustered.communities$vcount,
              "   COMMUNITIES#:", max(x$clustered.communities$membership),
              "   EDGES#:", ecount(x$clustered.network))

  
  i <- paste( attributes(x)$names, collapse=", ")
  #cat(l1, l2,"\n");cat(strwrap(paste(" + attr:",i,"\n"),width=60, exdent=8), sep="\n")
  cat(strwrap(paste(" + attr:",i,"\n"),width=60, exdent=8), sep="\n")
  cat(l1,"\n")
  print.igraph(x$raw.network)
  cat(l2,"\n")
  print.igraph(x$clustered.network)

}
