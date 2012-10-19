"print.pdb" <- function(x, ...) {
  
  ## This is adapted from the old "atom.select()" and
  ## "pdb.summary()" functions.

  eleno <- as.numeric(x$atom[,"eleno"])
  elety <- table(x$atom[,"elety"])
  resno <- rle2(as.numeric(x$atom[,"resno"]))
  resid <- table(x$atom[resno$inds,"resid"])
  segid <- unique(x$atom[,"segid"])
  chain <- unique(x$atom[,"chain"])
  ##chain <- rle2(x$atom[,"chain"])
  s <- paste(aa321(x$atom[resno$inds,"resid"]),collapse="")

  cat("..| segid |..",sep="\n");
  cat("total #: ")
  cat(length(segid),"\n")
  cat("values : ")
  cat(segid,"\n\n\n")

  cat("..| chain |..",sep="\n");
  cat("total #: ")
  cat(length(chain),"\n")
  cat("values : ")
  cat(chain,"\n\n\n")

  cat("..| resno |..",sep="\n");
  cat("total #: ")
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

