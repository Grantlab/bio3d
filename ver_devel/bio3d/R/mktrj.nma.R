"mktrj.nma" <- function(x=NULL,    # nma data structure
                        mode=7,      # which mode to move along
                        mag=10,      # magnification factor
                        step=1.25,   # step size
                        file=NULL,   # output pdb file
                        ... ) {      # args for write.pdb

  ## make a trjactory of atomic displacments along a given mode
  if(!"nma" %in% class(x))
    stop("mktrj.nma: must supply 'nma' object, i.e. from 'nma'")

  if(is.null(file))
    file <- paste("mode_", mode, ".pdb", sep="")

  if(x$L[mode]<=0)
    stop("Mode with eigenvalue <=0 detected. Check 'mode' index.")

  nstep <- c(seq(step, to=mag, by=step))
  zcoor <- cbind(1) %*% nstep

  scor  <- function(x,u,m) { return(x*u+m) }
  plus  <- sapply(c(zcoor), scor, u=x$modes[,mode], m=x$xyz)
  minus <- sapply(c(-zcoor), scor, u=x$modes[,mode], m=x$xyz)

  coor  <- t(cbind(x$xyz,
                   plus, plus[,rev(1:ncol(plus))],
                   x$xyz,
                   minus, minus[,rev(1:ncol(minus))]))

  write.pdb(xyz=coor, file=file, ...)
  invisible(coor)
}
