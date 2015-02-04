"aa321" <-
function (aa) {

  # convert three-letters amino-acid code into
  # one-letter IUPAC code, for instance "ALA" into "A".

  aa1 <- c("-", ".", "X", aa.mass$aa1, "H", "R")
  aa3 <- c("---", "---","UNK", aa.mass$aa3, "DDE", "CIR")

    convert <- function(x) {
      if(is.na(x)) return(NA)
      if (all(x != aa3)) {
        warning(paste("Unknown 3-letters code for aminoacid:",x))
        return("X") # mask unk
      }
      else {
        return(aa1[which(x == aa3)])
      }
    }
  return(as.vector(unlist(sapply(aa, convert))))
}


".aa321.na" <-
function (aa) {

  # convert three-letters amino-acid code into
  # one-letter IUPAC code, for instance "ALA" into "A".

  aa1 <- c("-",".","X",
           "C","G","T","A",
           "C","G","U","A")
  aa3 <- c("---", "---","UNK",
           "DC", "DG", "DT", "DA",
           "C",   "G", "T",  "A")
  
    convert <- function(x) {
      if(is.na(x)) return(NA)
      if (all(x != aa3)) {
        warning(paste("Unknown 3-letters code for aminoacid:",x))
        return("X") # mask unk
      }
      else {
        return(aa1[which(x == aa3)])
      }
    }
  return(as.vector(unlist(sapply(aa, convert))))
}

