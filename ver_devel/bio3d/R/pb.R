.init.pb <- function(ncore, min=0, max=1) {
  if(ncore == 1) {

     return ( txtProgressBar(min=min, max=max, style=3) )

  } else if(ncore > 1) {

     mcparallel <- get("mcparallel", envir = getNamespace("parallel"))
     mccollect <- get("mccollect", envir = getNamespace("parallel"))

     fpb <- fifo(tempfile(), open = "w+b", blocking = T)

     # spawn a child process for message printing
     child <- mcparallel({
        pb <- txtProgressBar(min=min, max=max, style=3)
        progress <- 0.0
        while(progress < max && !isIncomplete(fpb)) {
           msg <- readBin(fpb, "double")
           progress <- progress + as.numeric(msg)
           setTxtProgressBar(pb, progress)
        }
        close(pb)
     } )

     names(fpb) <- child$pid
     return(fpb)
  }
}
.update.pb <- function(pb, step=1) {

  if(inherits(pb, "txtProgressBar")) {
     i <- getTxtProgressBar(pb)
     setTxtProgressBar(pb, i+step)
  }
  else {
     if(inherits(pb, "fifo"))
        writeBin(step, pb)
  }

}
.close.pb <- function(pb) {
   if(inherits(pb, "fifo")) {
      mccollect <- get("mccollect", envir = getNamespace("parallel"))
      mccollect(as.numeric(names(pb)))
#      mccollect(as.numeric(names(pb)))
   }
   close(pb)
}
