summary.pdb <- function(object, ...) {

  ## summary.cna(pdb)

  if( !"pdb" %in% class(object) ) {
    stop("Input should be a pdb object, as obtained from 'read.pdb()'")
  }

  natom <- length(object$xyz)/3
  nprot <-length(atom.select(object, "protein", verbose=FALSE)$atom)
  nres <- sum(object$calpha)

  nhet <- nrow(object$het)
  if(is.null(nhet)) {
    nhet <- 0
    hetres <- "none"
    hetnres <- 0
  } else { 
    hetres <- paste( unique(object$het[,"resid"]), collapse=" ")
    hetnres <- length( unique(object$het[,"resno"]) )
  }
  
  not.prot.inds <- atom.select(object, "notprotein",verbose=FALSE)$atom
  nother.atom <- length( not.prot.inds )
  not.prot.res <- object$atom[not.prot.inds, "resid"]
  if(length(not.prot.res) == 0)
    not.prot.res <- "none"

  cat("\n Call:  ", paste(deparse(object$call), sep = "\n", collapse = "\n"),
        "\n", sep = "")

  s <- paste("\n  Atom Count:", natom+nhet,
             "\n\n   Total ATOMs#:", natom,
             "\n     Protein ATOMs#:", nprot,
             "  ( Calpha ATOMs#:", nres,")",
             "\n     Non-protein ATOMs#:", nother.atom,
             "  ( residues:", not.prot.res,")",
             "\n\n   Total HETATOMs:", nhet,
             "\n     Residues HETATOMs#:", hetnres,
             "  ( residues:", hetres,")\n")
              
  cat(s, "\n")
  i <- paste( attributes(object)$names, collapse=", ")
  cat(strwrap(paste(" + attr:",i,"\n"),width=45, exdent=8), sep="\n")
}
