
read.crd.amber2 <- function(file) {
  cl <- match.call()

  if(missing(file)) {
    stop("read.pdb: please specify a PDB 'file' for reading")
  }
  
  ##- Check if file exists locally or on-line
  if(!file.exists(file)) {
    stop("No input PDB file found: check filename")
  }
  
  ## parse CRD file with cpp function
  crd <- .read_crd(file)
  if(!is.null(crd$error))
    stop(paste("Could not read", file))
  else
    class(crd) <- c("amber", "crd")

  crd$xyz <- as.xyz(crd$xyz)
  crd$call <- cl
  
  ## finished
  return(crd)
}
