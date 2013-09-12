"aa321" <-
function (aa) {
  
  # convert three-letters amino-acid code into
  # one-letter IUPAC code, for instance "ALA" into "A".

  aa1 <- c("-",".","X",
           "A","C","D","E","F","G",
           "H","I","K","L","M","N","P","Q",
           "R","S","T","V","W","Y",
           "S","T","K","M","D","C","C","C","C","C","C","C","C"
           "H","H","H","H","H","H","H",
           "M", "D")
  aa3 <- c("---", "---","UNK",
           "ALA", "CYS", "ASP", "GLU", "PHE", "GLY",
           "HIS", "ILE", "LYS", "LEU", "MET", "ASN", "PRO", "GLN",
           "ARG", "SER", "THR", "VAL", "TRP", "TYR",
           "SEP", "TPO", "MLY", "MSE", "IAS", "ABA","CSO","CSD","CYM","CME","CSX","CMT","CYX"
           "HIE", "HIP", "HID", "HSD", "HSE", "HSP","DDE",
           "MHO", "ASX")
    convert <- function(x) {
      if(is.na(x)) return(NA)
      if (all(x != aa3)) {
        warning(paste("Unknown 3-letters code for aminoacid:",x))
        print(x)
        #return(NA)
        return("X") # mask unk
      }
      else {
        return(aa1[which(x == aa3)])
      }
    }
  return(as.vector(unlist(sapply(aa, convert))))
}

