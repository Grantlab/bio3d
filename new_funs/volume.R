
## Functions for outputing grid mapped volume data 
## to NetCDF format files suitable for reading with
## chimeria from UCSF

library(ncdf)

nt <- function(name, #units,
               dim, #missval,
               longname = name, prec = "single") {

  ## This is a modifed version of var.def.ncdf from the ncdf package

    oops <- require(ncdf)
    if (!oops) 
        stop("Please install the ncdf package from CRAN")
  
    if (!is.character(name)) {
        stop("Passed a var name that is NOT a string of characters!")
    }
#    if (storage.mode(missval) == "character") 
#        prec <- "char"
    var <- list()
    var$name <- name
    #var$units <- units
    #var$missval <- missval
    var$longname <- longname
    var$id <- -1
    attr(var, "class") <- "var.ncdf"
    if (prec == "float") 
        prec <- "single"
    if ((prec != "short") && (prec != "single") && (prec != "double") && 
        (prec != "integer") && (prec != "char") && (prec != "byte")) 
        stop(paste("var.def.ncdf: error: unknown precision specified:", 
            prec, ". Known values: short single double integer char byte"))
    var$prec <- prec
    if (is.character(class(dim)) && (class(dim) == "dim.ncdf")) 
        dim <- list(dim)
    var$dim <- dim
    var$ndims <- length(var$dim)
    for (i in 1:var$ndims) {
        if (class(var$dim[[i]]) != "dim.ncdf") {
            print(var)
            stop("Error, passed variable has a dim that is NOT of class dim.ncdf!")
        }
    }
    varunlimited <- FALSE
    if (var$ndims != 0) {
        for (i in 1:var$ndims) {
            if (var$dim[[i]]$unlim) 
                varunlimited <- TRUE
        }
    }
    var$unlim <- varunlimited
    return(var)
}


write.chimera.ncdf <- function(occupancy,
                               file="occupancy.nc") {

  ## Restrict to zero centered cube
  ## (i.e. box center at origin)
  origin=c(-150, -150, -150)
  full.dim <- dim(occupancy)
  cube <- all(full.dim == full.dim[1])

  if(!cube)
    stop("Currently we only deal with cubic grids")

  trim=FALSE
  
  if (trim) {
    ## trim occupancy matrix to remove zero full fields
    ind <- which(occupancy > 0, arr.ind=TRUE)

    center <- full.dim[1]/2
    x.unk <- unique(ind[,1])
    y.unk <- unique(ind[,2])
    z.unk <- unique(ind[,3])

    ## but insure we still have a cube
    x.bound <- max(abs( range(x.unk) - center))
    y.bound <- max(abs( range(y.unk) - center))
    z.bound <- max(abs( range(z.unk) - center))
    
    x.ind <- c(center-x.bound):c(x.bound+center)
    y.ind <- c(center-y.bound):c(y.bound+center)
    z.ind <- c(center-z.bound):c(z.bound+center)
    
    ##x.ind <- c( min(x.unk):max(x.unk))
    ##y.ind <- c( min(y.unk):max(y.unk))
    ##z.ind <- c( min(z.unk):max(z.unk))
    
    d.trim <- occupancy[x.ind, y.ind, z.ind]
    center <- dim(d.trim)/2 ## !! NOT FINISHED HERE YET !!
      
  } else {
    d.trim <- occupancy
  }

  ## Add one addational blank element per dim
  d.dim <- dim(d.trim)+1
  d.data <- array(0, dim = d.dim)
  d.data[-d.dim[1],-1,-1] <- d.trim


  ## Define xyz netcdf coordinate variables -- 1:dim
  dim1 = dim.def.ncdf( name="x", units="",
    1:d.dim[1], create_dimvar=FALSE)
  dim2 = dim.def.ncdf( name="y", units="",
    1:d.dim[2], create_dimvar=FALSE)
  dim3 = dim.def.ncdf( name="z", units="",
    1:d.dim[3], create_dimvar=FALSE)


  ## define the EMPTY (occupancy) data netcdf variable
  ##varz = var.def.ncdf(name="data", units="angstrom",
  ##             dim=list(dim1, dim2, dim3),  
  ##             missval = 1e+30, longname="data")
  ##source("new.def.R")
  
  varz = nt(name="data", dim=list(dim1, dim2, dim3))

  ## associate the netcdf variable with a netcdf file   
  nc.ex = create.ncdf( file, varz )

  ## add  global attributes
  att.put.ncdf( nc.ex, 0, "xyz_origin", origin )
  att.put.ncdf( nc.ex, 0, "xyz_step", c(1,1,1))

  ## put the variable into the file, and close
  put.var.ncdf(nc.ex, varz, d.data)

  ## add addational data attributes
  ##att.put.ncdf( nc.ex, "data", "rgba", c(0.7, 0.7, 0.7, 1.0) )
  ##att.put.ncdf( nc.ex, "data", "component_number", 1 )

  close.ncdf(nc.ex)

}

read.chimera.ncdf <- function(file) {
  ## file="occupancy.nc"
  nc <- open.ncdf(file, readunlim = FALSE)

  print(nc)
  ## summary(nc)

  ## x = get.var.ncdf( nc, "x")          # coordinate variable
  ## y = get.var.ncdf( nc, "y")          # coordinate variable
  ## z = get.var.ncdf( nc, "z")
  data = get.var.ncdf( nc, "data")   # variable

  origin <- att.get.ncdf( nc, 0, "xyz_origin")
  step  <- att.get.ncdf( nc, 0, "xyz_step")

  cat(paste("Global vars origin =", origin),sep="\n")
  cat(paste("Global vars step =", step),sep="\n")
  
  ##data.rgba <- att.get.ncdf( nc, "data", "rgba")
  ##data.comp <- att.get.ncdf( nc, "data", "component_number")

  close.ncdf(nc)
  return(data)
}
