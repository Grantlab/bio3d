print.prmtop <- function(x, printseq=TRUE, ...) {
  if(!is.null(x$SOLVENT_POINTER))
    sbox <- TRUE
  else
    sbox <- FALSE
  
  cn <- class(x)
  natom <-  x$POINTERS[1]
  ##nca <- length(which(x$ATOM_NAME=="CA" & x$AMBER_ATOM_TYPE=="CX"))
  
  if(sbox) {
    nres <- x$POINTERS[12]
    nres.solute <- x$SOLVENT_POINTER[1]

    nmol <- x$SOLVENT_POINTER[2]
    nmol.solute <- x$SOLVENT_POINTER[3]-1

    natom.per.mol <- x$ATOMS_PER_MOLECULE[1:nmol.solute]
    box.dim <- x$BOX_DIMENSIONS
  }
  else {
    nres <- x$POINTERS[12]
    nres.solute <- nres

    nmol <- length(which(x$ATOM_NAME=="OXT"))
    nmol.solute <- nmol
  }

  
  cat("\n Call:\n  ", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  
  cat(" Class:\n  ", paste(cn, collapse=", "),
      "\n\n", sep = "")

  cat(" System information:", "\n")
  cat("      Atom count:  ", natom, "\n", sep="")
  
  cat("      Residues:  ", nres.solute, " / ", nres, "\n", sep="")
  cat("      Molecules:  ", nmol.solute, " / ", nmol, "\n", sep="")

  if(sbox)
    cat("      Box dimensions:  ", paste(round(box.dim,2), collapse=" x "), "\n", sep="")


  if(printseq) {
    aa <- aa321(x$RESIDUE_LABEL)
    if(nres > 225) {
      ## Trim long sequences before output
      aa <- c(aa[1:225], "...<cut>...", aa[(nres-3):nres])
    }
    aa <- paste("     ",  gsub(" ","", 
            strwrap( paste(aa,collapse=" "), 
            width=120, exdent=0) ), collapse="\n")
    cat("\n")
    cat(" Sequence:\n", aa, "\n", sep="")
  }

  cat("\n")

  #i <- paste( attributes(x)$names, collapse=", ")
  #cat(strwrap(paste(" + attr:",i,"\n"),width=45, exdent=8), sep="\n")

  invisible( c(natom=natom, nres=nres, nres.solute=nres.solute,
               nmol=nmol, nmol.solute=nmol.solute) )
}
