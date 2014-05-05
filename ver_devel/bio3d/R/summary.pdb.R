summary.pdb <- function(object, printseq=FALSE, ...) {

  ## Print a summary of basic PDB object features

  if( !is.pdb(object) ) {
    stop("Input should be a pdb object, as obtained from 'read.pdb()'")
  }

  ## Multi-model check and total atom count
  nmodel <- nrow(object$xyz)
  if( is.null(nmodel) ) {
    ntotal <- length(object$xyz)/3
    nmodel = 1
  } else {
    ntotal <- length(object$xyz[1,])/3
  }

  nxyz <- length(object$xyz)
  nres <- sum(object$calpha)
  chains <- unique(object$atom[,"chain"])

  nprot <-length(atom.select(object, "protein", verbose=FALSE)$atom)
  other.inds <- atom.select(object, "notprotein", verbose=FALSE)$atom
  het <- object$atom[other.inds,]
  nhet.atom <- nrow(het)

  if(is.null(nhet.atom)) {
    nhet.atom <- 0
    nhet.res <- 0
    hetres <- "none"
  } else { 
  	hetres.resno <- apply(het[,c("chain","resno","resid")], 1, paste, collapse=".")
  	nhet.res <- length(unique(hetres.resno))
	hetres.nres <- table(het[,c("resid")][!duplicated(hetres.resno)])
  	hetres <- paste( paste0( names(hetres.nres), " (",hetres.nres, ")"), collapse=", ")
  }

  if((nprot+nhet.atom) != ntotal)
    warning("nPROTEIN + nNON-PROTEIN != nTotal")

 
  cat("\n Call:  ", paste(deparse(object$call), sep = "\n", collapse = "\n"),
        "\n", sep = "")

  s <- paste0("\n   Total Models#: ", nmodel, 
  			  "\n     Total Atoms#: ", ntotal, ",  XYZs#: ", nxyz, 
  			  "  Chains#: ", length(chains),
  			  "  (values: ", paste(chains, collapse=" "),")",

             "\n\n     Protein Atoms#: ", nprot,
             "  (residues/Calpha atoms#: ", nres,")",

             "\n\n     Non-protein Atoms#: ", nhet.atom,
             "  (residues: ", nhet.res, ")",
             "\n     Non-protein resid values: [", hetres,"]",
             "\n\n")
              
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

  invisible( c(nmodel=nmodel, natom=ntotal, nxyz=nxyz, nchains=length(chains),  
               nprot=nprot, nprot.res=nres, nother=nhet.atom, nother.res=nhet.res) )

}
