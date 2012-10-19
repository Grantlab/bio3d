`is.gap` <-
function(x, gap.char=c("-",".")) {
  return( as.logical( is.na(x) + (x %in% gap.char) ) )
}

