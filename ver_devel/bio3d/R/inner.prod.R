"inner.prod" <- function(a,b,mass=NULL) {
  if(!is.null(mass)) {
    if (nrow(as.matrix(a)) != (length(mass)*3)) 
      stop("unequal vector lengths")
  }
  
  if(is.null(mass))
    mass <- 1
  else
    mass <- rep(mass,each=3)

  if(class(a)=='matrix')
    return(colSums((a*b)*mass^2))
  else
    return(sum(a*b*mass^2))
}
