`read.ncdf` <-
function(trjfile, headonly = FALSE, verbose = TRUE){

  ##- Load ncdf package
  oops <- require(ncdf)
  if(!oops)
    stop("Please install the ncdf package from CRAN")
  
  nc <- open.ncdf(trjfile, readunlim=FALSE)

  conv <- att.get.ncdf( nc, varid=0, "Conventions")$value
  if(conv!="AMBER")
    warning("File conventions is not set to AMBER")

  if(verbose) {
    print(paste("Reading file", nc$filename))
    print(paste("Produced by program:",
      att.get.ncdf( nc, varid=0, "program")$value))
    print(paste("File conventions",conv, "version",
      att.get.ncdf( nc, varid=0, "ConventionVersion")$value))
    print(paste("Frames:",nc$dim$frame$len))
    print(paste("Atoms:", nc$dim$atom$len))
  }
  if(headonly) {
    ## Report only header information
    return(list("file"=nc$filename,
                "conv"=conv,
                "frames"=nc$dim$frame$len,
                "atoms"=nc$dim$atom$len))
  ##time  <- get.var.ncdf(nc,"time")
  }
  coords <- get.var.ncdf(nc,"coordinates")
  close.ncdf(nc)
  return(matrix(coords, ncol=(dim(coords)[2]*3),byrow=TRUE))
}

