`write.ncdf` <-
function(x, trjfile="R.ncdf"){

  ##- Load ncdf package
  oops <- require(ncdf)
  if(!oops)
    stop("Please install the ncdf package from CRAN")

  ## Error checking
  if(!is.matrix(x))
    stop("input x should be a natom by nframe numeric matrix of coordinates")

  nframe <- nrow(x)
  natom <- ncol(x)/3

  ## Define dimensions
  frame <- dim.def.ncdf(name="frame", units="", vals=c(1:nframe),
                        unlim=TRUE, create_dimvar=FALSE)

  spatial <- dim.def.ncdf(name="spatial", units="", vals=1:3, #c(1:3),#"xyz",
                          unlim=FALSE, create_dimvar=TRUE)

  atom <- dim.def.ncdf(name="atom", units="", vals=c(1:natom),
                       unlim=FALSE, create_dimvar=FALSE)

##  label <- dim.def.ncdf(name="label", units="", vals=1:5, ##???
##                        unlim=FALSE, create_dimvar=FALSE)

##  cells <- dim.def.ncdf(name="cell_spatial", units="", vals=1:3,# "abc",
##                        unlim=FALSE, create_dimvar=TRUE)

##  cella <- dim.def.ncdf(name="cell_angular", units="", vals=1:3,
##                        #vals=c("alpha", "beta", "gamma"),
##                        unlim=FALSE, create_dimvar=TRUE)

  ## Define variables
  time <- var.def.ncdf(name="time", units="picosecond", dim=frame,
                       missval=1e+30, prec="single") #"single" float
  coor <- var.def.ncdf(name="coordinates", units="angstrom", missval=1e+30,
                       dim=list(spatial,atom,frame), prec="single")#"single" float
##  cell.len <- var.def.ncdf(name="cell_lengths", units="angstrom", missval=1e+30,
##                           dim=list(cells,frame),  prec="double")
##  cell.ang <- var.def.ncdf(name="cell_angles", units="degrees", missval=1e+30,
##                           dim=list(cella,frame),  prec="double")

  ## Create the file
  ncw <- create.ncdf( trjfile, list(time, coor))#, cell.len, cell.ang) )

  ## Write data to file
  put.var.ncdf(ncw, time, c(1:nframe), start=1, count=nframe)
  put.var.ncdf( ncw, coor, t(x), start=c(1,1,1), count=c(3,natom,nframe))

  ## Define Required Attributes
  att.put.ncdf(ncw, varid=0, attname="Conventions", attval="AMBER")
  att.put.ncdf(ncw, varid=0, attname="ConventionVersion", attval="1.0")
  att.put.ncdf(ncw, varid=0, attname="program",attval="bio3d")
  att.put.ncdf(ncw, varid=0, attname="programVersion", attval="1.0-5")
  close.ncdf(ncw)
}

