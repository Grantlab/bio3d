library(bio3d)

## Description
## --------------------------
## Column 1: 3 letter aa code
## Column 2: 1 letter aa code
## Column 3: Formula collected from HICUP. 
##           e.g. http://xray.bmc.uu.se/hicup/CSD/
## Column 4: formula in peptide
## Column 5: Name of compound / amino acid

## Howto update 'aa.table':
## - Add new residue names at the bottom of 'rawtable'
## - source this file
## - move aa.table.rda to ../bio3d/data/aa.table.rda
## - commit changes

## Deviations from HICUP:
##  LYN, 

rawtable <- c(
  
  "ALA", "A", "C3 H7 N O2",     "C3 H5 N O1",    "Alanine", 
  "ARG", "R", "C6 H15 N4 O2",   "C6 H13 N4 O1",  "Arginine", 
  "ASN", "N", "C4 H8 N2 O3",    "C4 H6 N2 O2",   "Asparagine", 
  "ASP", "D", "C4 H6 N O4",     "C4 H4 N O3",    "Aspartic Acid", 
  "CYS", "C", "C3 H7 N O2 S",   "C3 H5 N O1 S",  "Cystein", 
  "GLN", "Q", "C5 H10 N2 O3",   "C4 H9 N2 O2",   "Glutamine", 
  "GLU", "E", "C5 H8 N O4",     "C5 H6 N O3",    "Glutamic Acid", 
  "GLY", "G", "C2 H5 N O2",     "C2 H3 N O1",    "Glycine", 
  "HIS", "H", "C6 H9 N3 O2",    "C6 H7 N3 O1",   "Histidine", 
  "ILE", "I", "C6 H13 N O2",    "C6 H11 N O1",   "Isoleucine", 
  "LEU", "L", "C6 H13 N O2",    "C6 H11 N O1",   "Leucine", 
  "LYS", "K", "C6 H15 N2 O2",   "C6 H13 N2 O1",  "Lysine", 
  "MET", "M", "C5 H11 N O2 S",  "C5 H9 N O1 S",  "Methionine", 
  "PHE", "F", "C9 H11 N O2",    "C9 H9 N O1",    "Phenylalanine", 
  "PRO", "P", "C5 H9 N O2",     "C5 H7 N O1",    "Proline", 
  "SER", "S", "C3 H7 N O3",     "C3 H5 N O2",    "Serine", 
  "THR", "T", "C4 H9 N O3",     "C4 H7 N O2",    "Threonine", 
  "TRP", "W", "C11 H12 N2 O2",  "C11 H10 N2 O1", "Tryptophan", 
  "TYR", "Y", "C9 H11 N O3",    "C9 H9 N O2",    "Tyrosine", 
  "VAL", "V", "C5 H11 N O2",    "C5 H9 N O1",    "Valine", 
  "HID", "H", "C6 H9 N3 O2",    "C6 H7 N3 O1",   "Histidine", 
  "HIE", "H", "C6 H9 N3 O2",    "C6 H7 N3 O1",   "Histidine", 
  "HIP", "H", "C6 H10 N3 O2",   "C6 H8 N3 O1",   "Histidine Positive", 
  "HSD", "H", "C6 H9 N3 O2",    "C6 H7 N3 O1",   "Histidine", 
  "HSE", "H", "C6 H9 N3 O2",    "C6 H7 N3 O1",   "Histidine", 
  "HSP", "H", "C6 H10 N3 O2",   "C6 H8 N3 O1",   "Histidine Positive", 
  "CYX", "C", "C3 H6 N O2 S",   "C3 H4 N O1 S",  "Cystein SSbond", 
  "ASH", "D", "C4 H7 N O4",     "C4 H5 N O3",    "Aspartic acid Neutral", 
  "GLH", "E", "C5 H9 N O4",     "C5 H7 N O3",    "Glutatmic acid Neutral", 
  "LYN", "K", "C6 H14 N2 O2",   "C6 H13 N2 O1",  "Lysine Neutral",
  "CYM", "C", "C3 H6 N O2 S",   "C3 H4 N O1 S",  "Cystein Negative", 
  "TPO", "T", "C4 H10 N O6 P",  "C4 H8 N O5 P",  "phosphothreonine", 
  "SEP", "S", "C3 H8 N O6 P",   "C3 H6 N O5 P",  "phosphoserine", 
  "PTR", "Y", "C9 H12 N O6 P",  "C9 H10 N O5 P", "o-phosphotyrosine", 
  "MLY", "K", "C8 H18 N2 O2",   "C8 H16 N2 O1",  "n-dimethyl-lysine", 
  "MSE", "M", "C5 H11 N O2 SE", "C5 H9 N O1 SE", "selenomethionine", 
  "IAS", "D", "C4 H7 N O4",     "C4 H5 N O3",    "beta-aspartyl", 
  "ABA", "X", "C4 H9 N1 O2",    "C4 H7 N1 O1",   "alpha-aminobutyric acid", 
  "CSO", "C", "C3 H7 N O3 S",   "C3 H5 N O2 S",  "s-hydroxycysteine", 
  "CSD", "C", "C3 H7 N O4 S",   "C3 H4 N O3 S",  "s-cysteinesulfinic acid", 
  "CME", "C", "C5 H11 N O3 S2", "C5 H9 N O2 S2", "s,s-(2-hydroxyethyl)thiocysteine", 
  "CSX", "C", "C3 H7 N O3 S",   "C3 H5 N O2 S",  "s-oxy cysteine", 
  "CMT", "C", "C4 H9 N O2 S",   "C4 H5 N O1 S",  "o-methylcysteine", 
  "MHO", "M", "C5 H11 N O3 S",  "C5 H9 N O2 S",  "s-oxymethionine", 
  "PFF", "F", "C9 H10 F N O2",  "C9 H8 F N O1",  "4-fluoro-l-phenylalanine", 
  "KCX", "K", "C7 H14 N2 O4",   "C7 H12 N2 O3",  "lysine nz-carboxylic acid", 
  "CSW", "C", "C3 H7 N O4 S",   "C3 H5 N O3 S",  "cysteine-s-dioxide", 
  "OCS", "C", "C3 H7 N O5 S",   "C3 H5 N O4 S",  "cysteinesulfonic acid", 
  "DDE", "H", "C13 H24 N5 O3",  "C13 H22 N5 O2", "diphthamide", 
  "CIR", "R", "C6 H13 N3 O3",   "C6 H11 N3 O2",  "citrulline"
  
  )

 
## convert to data frame
data <- as.data.frame(matrix(rawtable, ncol=5, byrow=TRUE), stringsAsFactors=FALSE)
colnames(data) <- c("aa3", "aa1", "formula_free", "formula", "name")

## calculate mass from formula
mass <- round(unlist(lapply(data$formula, formula2mass)), 3)

## add mass data
data <- cbind(data, mass)
data <- data[,c(1,2,6,4,5)]

## sort rows (standard amino acids first in table)
inds.a <- order(data$aa3[1:20])
inds.b <- order(data$aa3[21:nrow(data)])
inds <- c(inds.a, inds.b+20)

data <- data[inds, ]
rownames(data) <- data$aa3


## write matrix to file
##write.table(new.mat, quote=T, file="aa_table.mat.new", row.names=TRUE)

## convert to .rda
aa.table <- data
save(aa.table, file="aa.table.rda")



