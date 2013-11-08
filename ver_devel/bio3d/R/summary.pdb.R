summary.pdb <- function(object, printseq=FALSE, ...) {

  ## Print a summary of basic PDB object features

  if( !"pdb" %in% class(object) ) {
    stop("Input should be a pdb object, as obtained from 'read.pdb()'")
  }

  natom <- length(object$xyz)/3
  nprot <-length(atom.select(object, "protein", verbose=FALSE)$atom)
  nres <- sum(object$calpha)
  atom.chains <- unique(object$atom[,"chain"])

  nhet <- nrow(object$het)
  if(is.null(nhet)) {
    nhet <- 0
    hetres <- "none"
    hetnres <- 0
    het.chains <- 0
  } else { 
    hetres <- paste( unique(object$het[,"resid"]), collapse=" ")
    hetnres <- length( unique(object$het[,"resno"]) )
    het.chains <- unique(object$het[,"chain"])
  }
  
  ntotal <- natom+nhet
  not.prot.inds <- atom.select(object, "notprotein",verbose=FALSE)$atom
  nother.atom <- length( not.prot.inds )
  not.prot.res <- object$atom[not.prot.inds, "resid"]
  not.prot.res <- paste(unique(not.prot.res), collapse=" ")
 
  if(length(not.prot.res) == 0)
    not.prot.res <- "none"

  cat("\n Call:  ", paste(deparse(object$call), sep = "\n", collapse = "\n"),
        "\n", sep = "")

  s <- paste("\n  Atom Count:", ntotal,
             "\n\n   Total ATOMs#:", natom,
             "\n     Protein ATOMs#:", nprot,
             "  ( Calpha ATOMs#:", nres,")",
             "\n     Non-protein ATOMs#:", nother.atom,
             "  ( residues:", not.prot.res,")",
             "\n     Chains#:", length(atom.chains), 
             "  ( values:", paste(atom.chains, collapse=" "),")",

             "\n\n   Total HETATOMs:", nhet,
             "\n     Residues HETATOMs#:", hetnres,
             "  ( residues:", hetres,")",
             "\n     Chains#:", length(het.chains), 
             "  ( values:", paste(het.chains, collapse=" "),")\n\n")
              
  cat(s)

  if(printseq) {
    aa <- paste("     ",  gsub(" ","", 
            strwrap( paste(pdbseq(a),collapse=" "), 
            width=120, exdent=0) ), collapse="\n")
    cat("   Sequence:\n", aa, "\n\n", sep="")
  }

  i <- paste( attributes(object)$names, collapse=", ")
  cat(strwrap(paste(" + attr:",i,"\n"),width=45, exdent=8), sep="\n")

  invisible( c(ntotal=ntotal, ntotal=natom, nprot=nprot, 
               ca=nres, nother=nother.atom, ligs=not.prot.res,
               nhet=nhet, nhetres=hetnres, hetres=hetres) )
}
