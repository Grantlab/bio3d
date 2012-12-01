"print.pdb" <- function(x, ...) {
  
  ## This is adapted from the old "atom.select()" and
  ## "pdb.summary()" functions.

  ## Update to seperate water from sequence report
  y <- trim.pdb(x, atom.select(x, "water", verbose=F))
  p <- trim.pdb(x, atom.select(x, "notwater", verbose=F))

  ## Report on non-protein, possible ligand molecules
  i <- atom.select(p, "notprotein", verbose=F)
  no <- rle2(as.numeric(p$atom[i$atom,"resno"]))
  values <- (p$atom[i$atom,"resid"])[no$inds]

  ## Report on components of full PDB and non-water atoms
  eleno <- as.numeric(x$atom[,"eleno"])
  elety <- table(x$atom[,"elety"])
  resno <- rle2(as.numeric(p$atom[,"resno"]))
  resid <- table(x$atom[resno$inds,"resid"])
  segid <- unique(x$atom[,"segid"])
  chain <- unique(p$atom[,"chain"])

  ## Non water aa sequence
  s <- paste(suppressWarnings(aa321(p$atom[resno$inds,"resid"])),collapse="")

  ## Print output to terminal
  cat("..| segid |..",sep="\n");
  cat("total #: ")
  cat(length(segid),"\n")
  cat("values : ")
  cat(segid,"\n\n\n")

  cat("..| chain |..",sep="\n");
  cat("total protein #: ")
  cat(length(chain),"\n")
  cat("values : ")
  cat(chain,"\n\n\n")

  cat("..| resno |..",sep="\n");
  cat("total protein #: ")
  cat(length(resno$values),"\n")
  cat("in segments: \n\n")
  print(bounds(resno$values, pre.sort=F))
  cat("\n\n")

  cat("..| resid |..",sep="\n");
  cat("total #: ")
  cat(length(resid), "( diff types from",length(resno$values),")\n")
  cat("values : \n")
  print(resid)
  cat("\n\n")

  cat("..| eleno |..",sep="\n");
  cat("total # : ")
  cat(length(eleno),"\n")
  cat("in segments: \n\n")
  print(bounds(eleno, pre.sort=F))
  cat("\n\n")

  cat("..| elety |..",sep="\n");
  cat("total # : ")
  cat(length(elety),"( diff types from",length(eleno),")\n")
  cat("values : \n")
  print(elety)
  cat("\n\n")

  ## number of water molecules
  cat("..| water |..",sep="\n");
  cat("total # : ")
  cat( length( grep("O",y$atom[,"elety"]) ),"\n")
  cat("\n\n")

  ## number of non-protein molecules
  if(length(lig) > 0) {
    cat("..| other (not protein or water) |..",sep="\n");
    cat("total # : ")
    cat( length(unique(values)),"\n")
    ##cat("values : \n")
    print(table(values))
    cat("(These will be marked as 'X' in the sequence below.)\n\n\n")
  }
  
  cat("..| summary |..",sep="\n");
  cat("atom     #: ")
  cat(nrow(x$atom),"\n")
  cat("xyz      #: ")
  cat(length(x$xyz),"\n")
  cat("calpha   #: ")
  cat(sum(x$calpha),"\n")
  cat("sequence #: ")
  cat(s,"\n")


}

