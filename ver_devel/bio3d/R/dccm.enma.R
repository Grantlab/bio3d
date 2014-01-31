"dccm.enma" <- function(x, ...) {
  enma <- x
  if(!class(enma)=="enma")
    stop("input should be an 'enma' object as obtained from 'nma.pdbs'")
  
  ##if(is.null(enma$full.nma))
  ##  stop(paste("incompatible 'enma' object. \n", 
  ##             " run 'nma.pdbs' with 'full=TRUE'"))

  if(any(is.na(enma$U.subspace)))
    arr <- FALSE
  else
    arr <- TRUE
  
  dims <- dim(enma$U.subspace)

  if(arr)
    all.dccm <- array(0, dim=c(dims[1L]/3, dims[1L]/3, dims[3L]))
  else
    all.dccm <- list()

  ## Initialize progress bar
  pb <- txtProgressBar(min=1, max=dims[3L], style=3)
  
  for( i in 1:dims[3L] ) {

    if(is.null(enma$full.nma)) {
      dummy.nma <- list(U=enma$U.subspace[,,i],
                        L=enma$L[i,],
                        force.constants=enma$L[i,],
                        triv.modes=0,
                        natoms=nrow(enma$U.subspace[,,i])/3)
      class(dummy.nma) <- "nma"

      invisible(capture.output(
        cm.tmp <- dccm.nma(dummy.nma, ...) ))
    }
    else {
      invisible(capture.output(
        cm.tmp <- dccm.nma(enma$full.nma[[i]], ...) ))
    }

    if(arr)
      all.dccm[,,i] <- cm.tmp
    else
      all.dccm[[i]] <- cm.tmp

    setTxtProgressBar(pb, i)
  }
  close(pb)

  if(arr) {
    avg <- apply(all.dccm, 1:2, mean)
    class(avg) <- c("matrix", "dccm")
  }
  else
    avg <- NULL
  
  out <- list(all.dccm=all.dccm, avg.dccm=avg)
  return(out)
}
