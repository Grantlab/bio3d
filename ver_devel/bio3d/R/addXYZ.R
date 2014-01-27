addXYZ <- function(lwd = 2, labels= TRUE, cex = 2){
  
  oops <- require(rgl)
  if(!oops)
    stop("Please install the rgl package from CRAN")
  
  seg.id <- segments3d(
    rbind(
      c(0,0,0), c(5,0,0),
      c(0,0,0), c(0,5,0),
      c(0,0,0), c(0,0,5)
    ),
    lwd = lwd
  )
  seg.id <-data.frame(id = seg.id, type = "xyz.seg")

  lab.id <- NULL
  if(labels){
    lab.id <- text3d(
      rbind(
        c(0,0,0) + c(5.0+1.0,0.0    ,0.0    ),
        c(0,0,0) + c(0.0    ,5.0+1.0,0.0    ),
        c(0,0,0) + c(0.0    ,0.0    ,5.0+1.0)
      ),
      texts=c("x","y","z"),
      cex = cex
    )
    lab.id <- data.frame(id = lab.id, type = "xyz.lab")
  }
  ids <- rbind(seg.id, lab.id)
  
  invisible(ids)
}