`plot.hmmer` <-
function(x, ...) {
  allowed <- c("phmmer", "hmmsearch", "jackhmmer")
  
  if(!any(inherits(x, allowed)))
    stop(paste("please provide the results of a hmmer search of type:",
               paste(allowed, collapse=", ")))

  plot.blast(x, ...)
}
