summary.pdb2 <- function(object, printseq=FALSE, ...) {

  ## Print a summary of basic PDB object features

  if( !"pdb" %in% class(object) ) {
    stop("Input should be a pdb object, as obtained from 'read.pdb()'")
  }

  ## multi-model check and total atom count
  nmodel <- nrow(object$xyz)
  if( is.null(mm) ) {
    ntotal <- length(object$xyz)/3
    nmodel = 1
  } else {
    ntotal <- length(object$xyz[1,])/3
  }

  nxyz <- length(object$xyz)
  nprot <-length(atom.select(object, "protein", verbose=FALSE)$atom)
  nres <- sum(object$calpha)
  atom.chains <- unique(object$atom[,"chain"])

  ## HETATM
  het <- subset(object$atom, object$atom$type=="HETATM")
  nhet <- nrow(het)
  if(is.null(nhet)) {
    nhet <- 0
    hetres <- "none"
    hetnres <- 0
    het.chains <- 0
  } else { 
    hetres <- paste( unique(het[,"resid"]), collapse=" ")
    hetnres <- length( unique(het[,"resno"]) )
    het.chains <- unique(het[,"chain"])
  }
  
  ## ATOM
  object$atom <- subset(object$atom, object$atom$type=="ATOM")
  natom <- nrow(object$atom)

  if((natom+nhet) != ntotal)
    warning("nATOMs + nHETATMs != nTotal")

  not.prot.inds <- atom.select(object, "notprotein",verbose=FALSE)$atom
  nother.atom <- length( not.prot.inds )
  not.prot.res <- object$atom[not.prot.inds, "resid"]
  not.prot.res <- paste(unique(not.prot.res), collapse=" ")
 
  if(length(not.prot.res) == 0)
    not.prot.res <- "none"

  cat("\n Call:  ", paste(deparse(object$call), sep = "\n", collapse = "\n"),
        "\n", sep = "")

  s <- paste("\n  Atom Count:", ntotal, "  ( Total ATOMs/HETATOMs )",
             "\n      Models:", nmodel,
             "\n   XYZ Count:", nxyz,
             "\n\n   Total ATOMs#:", natom,
             "\n     Protein ATOMs#:", nprot,
             "  ( Calpha atoms#:", nres,")",
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
    aa <- pdbseq(object)
    if(nres > 225) {
      ## Trim long sequences before output
      aa <- c(aa[1:225], "...<cut>...", aa[(nres-3):nres])
    }
    aa <- paste("     ",  gsub(" ","", 
            strwrap( paste(aa,collapse=" "), 
            width=120, exdent=0) ), collapse="\n")
    cat("   Sequence:\n", aa, "\n\n", sep="")
  }

  i <- paste( attributes(object)$names, collapse=", ")
  cat(strwrap(paste(" + attr:",i,"\n"),width=45, exdent=8), sep="\n")

  invisible( c(ntotal=ntotal, nmodel=nmodel, nxyz=nxyz, natom=natom,  
               nprot=nprot, ca=nres, nother=nother.atom, ligs=not.prot.res,
               nhet=nhet, nhetres=hetnres, hetres=hetres) )
}
