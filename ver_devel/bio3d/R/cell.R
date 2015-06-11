is.cell <- function(x)
  inherits(x, "cell")

cell <- function(...)
  UseMethod("cell")

cell.default <- function(abc = rep(NA_real_, 3), abg = rep(90.0, 3),
                         orient = "cbzy", sgroup = "P 1", ...){
  if(!is.numeric(abc) || length(abc) != 3)
    stop("'abc' must be a length 3 numeric vector")
  if(!is.numeric(abg) || length(abg) != 3)
    stop("'abg' must be a length 3 numeric vector")
  if(!is.character(orient) || length(orient) != 1)
    stop("'orient' must be a length 1 character vector")
  if(!is.character(sgroup) || length(sgroup) != 1)
    stop("'sgroup' must be a length 1 character vector")
  names(abc) <- c("a", "b", "c")
  names(abg) <- c("alpha", "beta", "gamma")
  obj <- list(abc = abc, abg = abg, orient = orient, sgroup = sgroup)
  class(obj) <- c("cell", class(obj))
  return(obj)
}

cell.pdb <- function(x, ...)
  return(x$cell)

"cell<-" <- function(x, value)
  UseMethod("cell<-")

"cell<-.pdb" <- function(x, value){
  if(!inherits(value, "cell"))
    stop("'value' must be an object of class 'cell'")
  x["cell"] <- list(value)
  return(x)
}

cellCoords <- function(x)
  UseMethod("cellCoords")

cellCoords.cell <- function(x) {
  Pabc <- switch(substr(x$orient, 1, 2),
                 "ab" = list(Pabc = c(1, 2, 3), S =  1),
                 "ac" = list(Pabc = c(1, 3, 2), S = -1),
                 "bc" = list(Pabc = c(2, 3, 1), S =  1),
                 "ba" = list(Pabc = c(2, 1, 3), S = -1),
                 "ca" = list(Pabc = c(3, 1, 2), S =  1),
                 "cb" = list(Pabc = c(3, 2, 1), S = -1) )
  
  Pxyz <- switch(substr(x$orient, 3, 4),
                 "xy" = list(Pxyz = c(1, 2, 3), S =  1),
                 "xz" = list(Pxyz = c(1, 3, 2), S = -1),
                 "yz" = list(Pxyz = c(3, 1, 2), S =  1),
                 "yx" = list(Pxyz = c(2, 1, 3), S = -1),
                 "zx" = list(Pxyz = c(2, 3, 1), S =  1),
                 "zy" = list(Pxyz = c(3, 2, 1), S = -1) )
  
  S <- c(1, 1, Pabc$S*Pxyz$S)
  Pabc <- Pabc$Pabc
  Pxyz <- Pxyz$Pxyz
  
  abc <- unlist(x$abc[Pabc])
  abg <- unlist(x$abg[Pabc])*pi/180
  
  M <- matrix(NA_real_, nrow = 3, ncol = 3)
  M[, 1] <- c(abc[1], 0.0, 0.0)
  M[, 2] <- c(abc[2]*cos(abg[3]), abc[2]*sin(abg[3]), 0.0)
  M[1, 3] <- abc[3]*(cos(abg[2]))
  M[2, 3] <- abc[3]*(cos(abg[1]) - cos(abg[2])*cos(abg[3]))/sin(abg[3])
  M[3, 3] <- cellVolume(x)/(abc[1]*abc[2]*sin(abg[3]))
  dimnames(M) <- list(c("x", "y", "z")[Pxyz], c("a", "b", "c")[Pabc])
  M <- S*M[c("x", "y", "z"), c("a", "b", "c")]
  
  return(M)
}

cellCoords.pdb <- function(x){
  x <- cell(x)
  if(!is.null(x)) {
    cellCoords(x)    
  } else {
    stop("'x' doesn't have a 'cell' component")    
  }
}

cellVolume <- function(x)
  UseMethod("cellVolume")

cellVolume.cell <- function(x){
  x$abc[1]*x$abc[2]*x$abc[3]*
    sqrt(1 - sum(cos(unlist(x$abg)*pi/180)^2) +
           2*cos(x$abg[1]*pi/180)*cos(x$abg[2]*pi/180)*cos(x$abg[3]*pi/180))
}

cellVolume.pdb <- function(x){
  x <- cell(x)
  if(!is.null(x)) {
    cellVolume(x)    
  } else {
    stop("'x' doesn't have a 'cell' component")    
  }
}

cellDensity <- function(x)
  UseMethod("cellDensity")

cellDensity.pdb <- function(x){
  if(!is.null(x)) {
    Na <- 6.02214129E23
    sum(aa2mass(x))/(cellVolume(x)*1E-24*Na)
#   cat(d, " g.cm-3\n", sep = "")
  } else {
    stop("'x' doesn't have a 'cell' component")    
  }
}
