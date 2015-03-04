
read.prmtop2 <- function(file) {
  cl <- match.call()

  if(missing(file)) {
    stop("read.pdb: please specify a PDB 'file' for reading")
  }
  
  ##- Check if file exists locally or on-line
  if(!file.exists(file)) {
    stop("No input PDB file found: check filename")
  }
  
  ## parse PRMTOP file with cpp function
  prmtop <- .read_prmtop(file)
  if(!is.null(prmtop$error))
    stop(paste("Could not read", file))
  else
    class(prmtop) <- c("amber", "prmtop")
  
  prmtop$call <- cl
  
  ## finished
  return(prmtop)
}
