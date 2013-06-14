`write.ncdf` <-
function(x, cell=NULL, trjfile="R.ncdf"){

  ##- Load ncdf package
  oops <- require(ncdf)
  if(!oops)
    stop("Please install the ncdf package from CRAN")

  ## Error checking
  if(!is.matrix(x))
    stop("input x should be a natom by nframe numeric matrix of coordinates")
  
  if(!is.null(cell)) {
     if(!is.matrix(cell))
       stop("input cell should be a 6 by nframe numeric matrix")
  }

  nframe <- nrow(x)
  natom <- ncol(x)/3

  ## Define dimensions
  frame <- dim.def.ncdf(name="frame", units="", vals=c(1:nframe),
                        unlim=TRUE, create_dimvar=FALSE)

  spatial <- dim.def.ncdf(name="spatial", units="", vals=1:3, #c(1:3),#"xyz",
                          unlim=FALSE, create_dimvar=TRUE)

  atom <- dim.def.ncdf(name="atom", units="", vals=c(1:natom),
                       unlim=FALSE, create_dimvar=FALSE)

  if(!is.null(cell)) {
     label <- dim.def.ncdf(name="label", units="", vals=1:5, 
             unlim=FALSE, create_dimvar=FALSE)
     cell_spatial <- dim.def.ncdf(name="cell_spatial", units="", 
            vals=1:3, unlim=FALSE, create_dimvar=TRUE)
     cell_angular <- dim.def.ncdf(name="cell_angular", units="", 
            vals=1:3, unlim=FALSE, create_dimvar=TRUE)
  }
     

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
  if(!is.null(cell)) {
    cell_lengths <- var.def.ncdf(name="cell_lengths", units="angstrom", 
                 missval=1e+30, dim=list(cell_spatial, frame), prec="double")
    cell_angles <- var.def.ncdf(name="cell_angles", units="degree",
                 missval=1e+30, dim=list(cell_angular, frame), prec="double")
  }
##  cell.len <- var.def.ncdf(name="cell_lengths", units="angstrom", missval=1e+30,
##                           dim=list(cells,frame),  prec="double")
##  cell.ang <- var.def.ncdf(name="cell_angles", units="degrees", missval=1e+30,
##                           dim=list(cella,frame),  prec="double")
  
  ## Create the file
  if(!is.null(cell)) {
     ncw <- create.ncdf( trjfile, list(time, coor, cell_lengths, cell_angles))
  } else {
     ncw <- create.ncdf( trjfile, list(time, coor))#, cell.len, cell.ang) )
  }

  ## Write data to file
  put.var.ncdf(ncw, time, c(1:nframe), start=1, count=nframe)
  put.var.ncdf( ncw, coor, t(x), start=c(1,1,1), count=c(3,natom,nframe))
  if(!is.null(cell)) {
     put.var.ncdf( ncw, cell_lengths, t(cell[,1:3]), start=c(1,1), count=c(3,nframe))
     put.var.ncdf( ncw, cell_angles, t(cell[,4:6]), start=c(1,1), count=c(3,nframe))
  }

  ## Define Required Attributes
  att.put.ncdf(ncw, varid=0, attname="Conventions", attval="AMBER")
  att.put.ncdf(ncw, varid=0, attname="ConventionVersion", attval="1.0")
  att.put.ncdf(ncw, varid=0, attname="program",attval="bio3d")
  att.put.ncdf(ncw, varid=0, attname="programVersion", attval="1.2")
  close.ncdf(ncw)
}

