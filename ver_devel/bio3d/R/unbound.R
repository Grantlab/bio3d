"unbound" <-
function(start,end) {
  if(length(start)!=length(end))
    stop("start and end must are not the same length")
  ex <- NULL
  for(i in 1:length(start)) {
    ex <- c(ex, start[i]:end[i])
  }
  return(ex)
}

