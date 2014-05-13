addPBCBox <- function(x, lwd = 2, col = "white"){

  oops <- require(rgl)
  if(!oops)
    stop("Please install the rgl package from CRAN")
  
  cell <- cell.coords(x)
  
  seg.id <- segments3d(
    rbind(
      c(0,0,0)                  , cell[,1],
      c(0,0,0)                  , cell[,2],
      c(0,0,0)                  , cell[,3],
      cell[,1]+cell[,2]         , cell[,1],
      cell[,1]+cell[,2]         , cell[,2],
      cell[,1]+cell[,2]         , cell[,1]+cell[,2]+cell[,3],
      cell[,1]+cell[,3]         , cell[,1],cell[,1]+cell[,3],
      cell[,1]+cell[,3]+cell[,2], cell[,1]+cell[,3],cell[,3],
      cell[,2]+cell[,3]         , cell[,2]+cell[,3]+cell[,1],
      cell[,2]+cell[,3]         , cell[,2],
      cell[,2]+cell[,3]         , cell[,3]
    ),
    col = col,
    lwd = lwd
  )
  seg.id <- data.frame(id = seg.id, type = "pbc.box")
  
  invisible(seg.id)
}