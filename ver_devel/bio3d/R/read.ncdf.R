`read.ncdf` <-
function(trjfile, headonly = FALSE, verbose = TRUE, time=FALSE,
         start = NULL, end = NULL){
  
  # Currently file open in SHARE mode is not supported by NCDF
  # Multicore support for reading single file is supressed 
  ncore = 1
  nseg.scale = 1
  # Parallelized by multicore package (Wed Apr 24 10:40:39 EDT 2013)
  if(ncore > 1) {
     oops <- require(multicore)
     if(!oops)
        stop("Please install the multicore package from CRAN")

     options(cores = ncore)

     # Issue of serialization problem
     # Maximal number of cells of a double-precision matrix
     # that each core can serialize: (2^31-1-61)/8
     R_NCELL_LIMIT_CORE = 2.68435448e8
     R_NCELL_LIMIT = ncore * R_NCELL_LIMIT_CORE

     if(nseg.scale < 1) {
        warning("nseg.scale should be 1 or a larger integer\n")
        nseg.scale=1
     }
  }

  ##- Load ncdf package
  oops <- require(ncdf)
  if(!oops)
    stop("Please install the ncdf package from CRAN")
  
  # open files 
  nc <- lapply(trjfile, function (fnm) {
#     if(verbose) print(paste("Opening file ", fnm, sep="")) 

     nc <- open.ncdf(fnm, readunlim=FALSE)
     conv <- att.get.ncdf( nc, varid=0, "Conventions")$value
     if(conv!="AMBER")
        warning(paste("File conventions is not set to AMBER",
                  fnm, sep=": "))
     return(nc)
  })
  
  # prepare frame start and end No for each file
  # for time range selection
  frange <- NULL
  if(!all(is.null(c(start, end))) && !time) {
     flen <- unlist(lapply(nc, function(n) return(n$dim$frame$len)))
     frange <- matrix(c(1, cumsum(flen[-length(nc)])+1, cumsum(flen)), 
                    nrow=length(nc))
  }
  
  # set range of frames for reading
  # skip files out of the selection range
  # read heads or coordinates
  retval <- lapply(1:length(nc), function (inc) {
     nc <- nc[[inc]]  
     if(!is.null(frange)) frange <- frange[inc,]
     ss = 1
     ee = nc$dim$frame$len
     if (!all(is.null(c(start, end)))) {
        #check frame No or time
        btime = frange[1]
        etime = frange[2]
        if(time) {
           btime <- get.var.ncdf(nc, "time", 1, 1)
           etime <- get.var.ncdf(nc, "time", nc$dim$frame$len, 1)
        }
        if((!is.null(start) && (etime < start)) || 
              (!is.null(end) && end >=0 && btime > end) ||
             (!is.null(start) && !is.null(end) && end >=0 && end < start) ) {
           if(verbose) print(paste("Skip file", nc$filename))
           close.ncdf(nc)
           return()
        } 
        if(!headonly) {
           timeall <- btime:etime
           if(time) timeall <- get.var.ncdf(nc, "time")
           ss <- if(is.null(start)) 1 else which((timeall - start) >=0 )[1]
           if(is.null(end) || end < 0 || end > etime) {
              ee = nc$dim$frame$len
           } else {
              ee <- which((timeall - end) <= 0) 
              ee <- ee[length(ee)]
           }
        }
     }
     tlen = ee - ss + 1
     if(verbose) {
       print(paste("Reading file", nc$filename))
       print(paste("Produced by program:",
         att.get.ncdf( nc, varid=0, "program")$value))
       conv <- att.get.ncdf( nc, varid=0, "Conventions")$value
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
     coords <- get.var.ncdf(nc, "coordinates", c(1, 1, ss), 
                             c(-1, -1, tlen))
     close.ncdf(nc)
     return(matrix(coords, ncol=(dim(coords)[2]*3),byrow=TRUE))
  } )

  if(headonly) {
     retval <- do.call("c", retval)
  } else {
     retval <- do.call(rbind, retval)
  }
  return(retval)
}
