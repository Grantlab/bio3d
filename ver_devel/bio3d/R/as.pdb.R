as.pdb <- function(...)
  UseMethod("as.pdb")


as.pdb.default <- function(xyz=NULL, type = NULL, resno = NULL,
                           resid = NULL, eleno = NULL,
                           elety = NULL, chain = NULL,
                           insert= NULL, alt = NULL,
                           o = NULL,  b = NULL,
                           segid = NULL, elesy = NULL,
                           charge = NULL, verbose=TRUE, ...) {
  cl <- match.call()

  ## determine number of atoms
  input <- list(xyz=xyz, eleno=eleno, resno=resno, resid=resid)
  nulls <- unlist(lapply(input, is.null))
  inds  <- which(!nulls)

  if(length(inds)==0)
    stop("insufficient arguments. provide 'xyz', 'eleno', 'resno', and/or 'resid'")

  ## if xyz is provided use it to determine natoms
  if (inds[1]==1) {
    xyz <- as.xyz(xyz)
    natoms <- ncol(xyz)/3
  }
  else {
    natoms <- length(input[[inds[1]]])
  }
  
  if(verbose) {
    cat("\n")
    cat(" Summary of PDB generation:\n")
    cat(paste(" .. number of atoms in PDB determined by '", names(input)[inds[1]], "'\n", sep=""))
  }

  ## set value of 'xyz'
  if(!is.null(xyz)) {
    if(!is.numeric(xyz) | !is.xyz(xyz))
      stop("'xyz' must be a numeric vector/matrix")
    xyz <- as.xyz(xyz)
    if((ncol(xyz)/3)!=natoms)
      stop("ncol(xyz)/3 != length(resno)")
  }
  else {
    xyz <- as.xyz(rep(NA, natoms*3))
  }

  ## generic function to set the values of remaining columns of PDB
  .setval <- function(values=NULL, typ=NULL, default=NULL, class="character", natoms=0, repfirst=FALSE) {
    if(!is.null(values)) {
      if(class=="character") fun=is.character
      if(class=="numeric") fun=is.numeric
      if(!fun(values))
        stop(paste("'", typ, "' must be a ", class, " vector", sep=""))
      
      if(length(values)==1 & repfirst)
        values <- rep(values, natoms)
      
      if(length(values)!=natoms)
        stop(paste("length(", typ, ") != natoms", sep=""))
    }
    else {
      values <- default
    }
    return(values)
  }

  type  <- .setval(type, typ="type", default=rep("ATOM", natoms), class="character", natoms=natoms, repfirst=TRUE)
  eleno  <- .setval(eleno, typ="eleno", default=seq(1, natoms), class="numeric", natoms=natoms, repfirst=FALSE)
  elety  <- .setval(elety, typ="elety", default=rep("CA", natoms), class="character", natoms=natoms, repfirst=TRUE)
  resno  <- .setval(resno, typ="resno", default=seq(1, natoms), class="numeric", natoms=natoms, repfirst=FALSE)
  chain  <- .setval(chain, typ="chain", default=rep("A", natoms), class="character", natoms=natoms, repfirst=TRUE)
  resid  <- .setval(resid, typ="resid", default=rep("ALA", natoms), class="character", natoms=natoms, repfirst=TRUE)
  elesy  <- .setval(elesy, typ="elesy", default=rep(NA, natoms), class="character", natoms=natoms, repfirst=TRUE)
  segid  <- .setval(segid, typ="segid", default=rep(NA, natoms), class="character", natoms=natoms, repfirst=TRUE)
  o  <- .setval(o, typ="o", default=rep(NA, natoms), class="numeric", natoms=natoms, repfirst=TRUE)
  b  <- .setval(b, typ="b", default=rep(NA, natoms), class="numeric", natoms=natoms, repfirst=TRUE)
  alt  <- .setval(alt, typ="alt", default=rep(NA, natoms), class="character", natoms=natoms, repfirst=FALSE)
  insert <- .setval(insert, typ="insert", default=rep(NA, natoms), class="character", natoms=natoms, repfirst=FALSE)
  charge  <- .setval(charge, typ="charge", default=rep(NA, natoms), class="numeric", natoms=natoms, repfirst=TRUE)

  ## make the data frame for the final PDB object
  atom        <- list()
  atom$type   <- type
  atom$eleno  <- eleno 
  atom$elety  <- elety
  atom$alt    <- alt
  atom$resid  <- resid
  atom$chain  <- chain
  atom$resno  <- resno
  atom$insert <- insert
  atom$x      <- xyz[1, seq(1, natoms*3, by=3)]
  atom$y      <- xyz[1, seq(2, natoms*3, by=3)]
  atom$z      <- xyz[1, seq(3, natoms*3, by=3)]
  atom$o      <- o
  atom$b      <- b
  atom$segid  <- segid
  atom$elesy  <- elesy
  atom$charge <- charge
  atom        <- data.frame(atom, stringsAsFactors=FALSE)

  out        <- list()
  out$atom   <- atom
  out$xyz    <- xyz
  class(out) <- "pdb"
  
  ca.inds    <- atom.select.pdb(out, "calpha", verbose=FALSE)
  out$calpha <- seq(1, natoms) %in% ca.inds$atom
  out$call <- cl

  if(verbose) {
    resid <- unique(paste(atom$chain, atom$resno, sep="-"))
    cat(paste(" .. number of atoms in PDB: ", natoms, "\n"))
    cat(paste(" .. number of calphas in PDB:", sum(out$calpha), "\n"))
    cat(paste(" .. number of residues in PDB:", length(resid), "\n"))
    cat("\n")
  }
  
  return(out)
}
