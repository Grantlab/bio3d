`read.ncdf` <-
function(trjfile, headonly = FALSE, verbose = TRUE, time=FALSE,
         first = NULL, last= NULL, stride = 1, cell = FALSE,
         at.sel = NULL){
  # Adding option 'at.sel' to select a subset of the structure using
  # an object of class 'select'
  
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
  
  # prepare frame first and last No for each file
  # for time range selection
  frange <- NULL
  if(!is.null(c(first, last)) && !time) {
     flen <- unlist(lapply(nc, function(n) return(n$dim$frame$len)))
     frange <- matrix(c(1, cumsum(flen[-length(nc)])+1, cumsum(flen)), 
                    nrow=length(nc))
  }
  
  # set range of frames for reading
  # skip files out of the selection range
  # read heads or coordinates
  retval <- lapply(seq_along(nc), function (inc) {
    first.atom <-  1
    count.atom <- -1
    if(!is.null(at.sel)) {
      if(!is.select(at.sel)) stop("'at.sel' must be an object of class 'select'. See 'atom.select'.")
      atom.ind <- xyz2atom(at.sel$xyz)
      first.atom <- min(atom.ind)
      count.atom <- diff(range(atom.ind)) + 1
    }
    
     nc <- nc[[inc]]  
     if(!is.null(frange)) frange <- frange[inc,]
     ss = 1
     ee = nc$dim$frame$len
     if (!is.null(c(first, last))) {

        #check frame No or time
        if(time) {
           btime <- get.var.ncdf(nc, "time", 1, 1)
           etime <- get.var.ncdf(nc, "time", nc$dim$frame$len, 1)
        } else {
           btime = frange[1]
           etime = frange[2]
        }

        if((!is.null(first) && (etime < first)) || 
              (!is.null(last) && last >=0 && btime > last) ||
             (!is.null(first) && !is.null(last) && last >=0 && last < first) ) {
           if(verbose) print(paste("Skip file", nc$filename))
           close.ncdf(nc)
           return()
        } 
        if(!headonly) {
           timeall <- btime:etime
           if(time) timeall <- get.var.ncdf(nc, "time")

           ss <- if(is.null(first)) 1 else which((timeall - first) >=0 )[1]
           if(is.null(last) || last < 0 || last > etime) {
              ee = nc$dim$frame$len
           } else {
              ee <- which((timeall - last) <= 0) 
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
     if(cell) {
        celldata <- get.var.ncdf(nc, "cell_lengths", c(1, ss), 
                             c(-1, tlen))
        celldata <- t( rbind(celldata, get.var.ncdf(nc, "cell_angles", 
                         c(1, ss), c(-1, tlen))) )
        if(time)
           rownames(celldata) <- get.var.ncdf(nc, "time", ss, tlen)
        close.ncdf(nc)
        return( celldata )
     }
     if(count.atom < 0) count.atom = nc$dim$atom$len

     # To solve 32-bit limitation problem with large trajectory file
     .check <- (3 * count.atom * tlen) / (2^31 - 1)
     if(.check > 1) {
        .nb <- floor(.check) + 1
        .nn <- floor(tlen / .nb)
        .ss <- seq(ss, ss + tlen - 1, .nn)
        .tlen <- rep(.nn, length(.ss) - 1)
        .tlen <- c(.tlen, tlen - sum(.tlen))
        coords <- sapply(1:length(.ss), function(i) 
            get.var.ncdf(nc, "coordinates", c(1, first.atom, .ss[i]), 
                          c(-1, count.atom, .tlen[i])) )
     } else {    
        coords <- get.var.ncdf(nc, "coordinates", c(1, first.atom, ss), 
                          c(-1, count.atom, tlen))
     }
     if(!is.null(at.sel)) coords <- coords[,atom.ind - first.atom + 1,]
     coords <- matrix( coords, ncol=(dim(coords)[2]*3), byrow=TRUE )
     if(time)
        rownames(coords) <- get.var.ncdf(nc, "time", ss, tlen)
     close.ncdf(nc)
     return( coords )
  } )

  if(headonly) {
     retval <- do.call("c", retval)
  } else {
     retval <- do.call(rbind, retval)
     # take every "stride" frame
     retval <- subset(retval, (1:nrow(retval)) %in% seq(1, nrow(retval), stride))
  }
  class(retval)="xyz"
  return(retval)
}
