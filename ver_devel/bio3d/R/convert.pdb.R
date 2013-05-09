"convert.pdb" <-
function(pdb, type,
         renumber=FALSE, first.resno=1, first.eleno=1,
         rm.h=TRUE, rm.wat=FALSE) {

  ## !!This function is very poor and needs to be re-writen!!
  
  options <- c("ori", "pdb","charmm","amber","gromacs")
  type <- match.arg(tolower(type), options)

  if(type != "ori") {
  cat(paste(" Converting to", type, "format"),"\n")

  ## residue and atom types from PDB
  restype <- unique(pdb$atom[,"resid"])
  #eletype <- unique(pdb$atom[,"elety"])

  ## Determine input PDB type based on restype and eletype
  ## NOT IMPLIMENTED YET!

  
  ## different his resids
  his <- matrix( c("HIS", "HSD", "HID","HISA",
                   "HIS", "HSE", "HIE","HISB",
                   "HIS", "HSP", "HIP","HISH"),
                nrow=3, byrow=TRUE,
                dimnames = list(c("d","e","b"),
                  c("pdb","charmm","amber","gromacs")) )
  
  standard <- c("MET", "PRO", "LEU", "VAL", "ALA",
                "ASP", "GLY", "ARG", "CYS", "THR",
                "GLU", "ILE", "LYS", "PHE", "GLN",
                "SER", "ASN", "TYR", "TRP", his)

  non.standard <- !(restype %in% standard)
  
  if(sum( non.standard ) > 0) {
    cat(paste(" Found ",sum( non.standard ),
              " non standard residues:"),
        restype[non.standard], sep="\n")
  }

  ## Convert HIS resid
  type.inds <-  (colnames(his) %in% type)
  conv.inds <- !(colnames(his) %in% c(type,"pdb"))

  his.d.ind <- (pdb$atom[,"resid"] %in% his["d", !type.inds ])
  his.e.ind <- (pdb$atom[,"resid"] %in% his["e", conv.inds ])
  his.b.ind <- (pdb$atom[,"resid"] %in% his["b", conv.inds ])
  
  pdb$atom[his.d.ind,"resid"] <- his["d", type.inds ]
  pdb$atom[his.e.ind,"resid"] <- his["e", type.inds ]
  pdb$atom[his.b.ind,"resid"] <- his["b", type.inds ]

  ## Convert ILE CD elety
  if (type=="charmm") {
    ile.ind <- colSums( rbind((pdb$atom[,"elety"] %in% "CD1"),
                              (pdb$atom[,"resid"] %in% "ILE"))) ==2
    pdb$atom[ile.ind,"elety"] <- "CD"

    pdb$atom[,"chain"]=NA # strip chain
  } else {
    ile.ind <- colSums( rbind((pdb$atom[,"elety"] %in% "CD"),
                            (pdb$atom[,"resid"] %in% "ILE"))) ==2
    pdb$atom[ile.ind,"elety"] <- "CD1"
  }
} ## END type != "ori"

  ## Strip hydrogen
  if(rm.h) {
    h.ind <- which(substr(pdb$atom[,"elety"],1,1) %in% "H")
    if(length(h.ind) > 0) {
      pdb$atom <- pdb$atom[-h.ind,]
      pdb$xyz <- pdb$xyz[-c((h.ind*3)-2, (h.ind*3)-1, (h.ind*3))]
      ##pdb$calpha <- pdb$calpha[-h.ind]
    }
    if(!is.null(pdb$het)) {
      h.ind <- which(substr(pdb$het[,"elety"],1,1) %in% "H")
      if(length(h.ind) > 0) {
        pdb$het <- pdb$het[-h.ind,]
      }
    }
  } else {
    ## convert hydrogen atom types
    if(type=="pdb") {
      pdb$atom[ pdb$atom[,"elety"]=="HN", "elety"] = "H"
      ###!!! DO MORE ATOM TYPE CONVERSION WHEN I GET TIME !!!###
    }
  }


  
  ## Strip water
  if(rm.wat) {
    wat <- c("H2O",  "OH2", "HOH", "HHO", "OHH", "SOL",
             "WAT", "TIP", "TIP2", "TIP3", "TIP4")
  
    wat.ind <- which(pdb$atom[,"resid"] %in% wat)
    if(length(wat.ind) > 0) {
      pdb$atom <- pdb$atom[-wat.ind,]
      pdb$xyz <- pdb$xyz[-c((wat.ind*3)-2, (wat.ind*3)-1, (wat.ind*3))]
      pdb$calpha <- pdb$calpha[-wat.ind]
    }
  }

  ## Renumber 
  if(renumber) {
    pdb$atom[,"eleno"] <- seq(first.eleno, length=nrow(pdb$atom))
    s.ind <- which(!duplicated(pdb$atom[,"chain"]))
    e.ind   <- c(s.ind[-1]-1, nrow(pdb$atom))

    ibase = 0
    for (i in 1:length(s.ind)) {
       nums <- as.numeric(pdb$atom[s.ind[i]:e.ind[i],"resno"])
       ##pdb$atom[,"resno"] <- nums - (nums[1] - first.resno)
   
       ## concetive residue numbers
       tbl <- table(nums)
       #new.nums <- first.resno:(first.resno+length(tbl))
       new.nums <- (first.resno+ibase):(first.resno+length(tbl)-1+ibase)
       pdb$atom[s.ind[i]:e.ind[i],"resno"] <- rep(new.nums, tbl)
       ibase = ibase + length(tbl)
    }
  }
  
  ## (split by chain or segid)
  #unique(pdb$atom[,"chain"])
  #unique(pdb$atom[,"segid"])

  return(pdb)
}

