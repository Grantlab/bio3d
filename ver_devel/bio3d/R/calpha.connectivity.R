calpha.connectivity <- function(...)
  UseMethod("calpha.connectivity")

calpha.connectivity.xyz <- function(x, d.cut = 4, ...){
  if(!inherits(x, "xyz"))
    stop("'x' must be an object of class 'xyz'")
  if(length(x) == 0)
    stop("'x' is empty")
  
  natom <- length(x)/3
  con.all <- connectivity(eleno.1 = 1:(natom-1), 
                          eleno.2 = 2: natom)
  
  a <- matrix(x[atom2xyz(con.all$eleno.1)], ncol=3, byrow=T)
  b <- matrix(x[atom2xyz(con.all$eleno.2)], ncol=3, byrow=T)
  d <- dist.xyz(a, b, all.pairs=FALSE)
  inds <- d < d.cut
  con <- con.all[inds,]
  con <- connectivity.default(con$eleno.1, con$eleno.2)
  return( con )
}

calpha.connectivity.pdb <- function(x, d.cut = 4, ...){
  if(!is.pdb(x))
    stop("'x' must be an object of class 'pdb'")
  
  x <- trim.pdb(x, atom.select(x, "calpha", verbose=FALSE))
  calpha.connectivity.xyz(x$xyz, d.cut = d.cut, ...)
}
