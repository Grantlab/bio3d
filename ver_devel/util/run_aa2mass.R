library(bio3d)

if(FALSE) {
w <- c( 71.080, 157.204, 114.108, 114.082, 103.150,
       128.134, 128.108,  57.054, 137.146, 113.158,
       113.158, 129.184, 131.202, 147.172,  97.116,
       87.080, 101.106, 186.210, 163.172,  99.132,
       
       137.146, 137.146, 138.154, 137.146, 137.146, 138.154,
       102.142, 115.09, 129.116, 128.176, 102.142,
       
       195.068, 181.042, 241.134, 156.228, 131.202, 101.106,
       85.106, 119.15,  134.142, 179.272, 119.150, 103.150,
       147.202, 165.164)
}

## Specify 3-letter AA code
aa <- c("ALA", "ARG", "ASN", "ASP", "CYS",
        "GLN", "GLU", "GLY", "HIS", "ILE",
        "LEU", "LYS", "MET", "PHE", "PRO",
        "SER", "THR", "TRP", "TYR", "VAL",
        
        "HID", "HIE", "HIP", "HSD", "HSE", "HSP",
        "CYX", "ASH", "GLH", "LYN", "CYM",
        
        "TPO", "SEP", "PTR", "MLY", "MSE", "IAS",
        "ABA", "CSO", "CSD", "CME", "CSX", "CMT",
        "MHO", "PFF", "KCX", "CSW", "OCS")

## Amino acid names
desc <- c("Alanine", "Arginine", "Asparagine", "Aspartic Acid", "Cystein",
          "Glutamine", "Glutamic Acid","Glycine","Histidine","Isoleucine",
          "Leucine","Lysine", "Methionine", "Phenylalanine","Proline",
          "Serine", "Threonine", "Tryptophan","Tyrosine","Valine",

          "Histidine", "Histidine", "Histidine Positive", "Histidine", "Histidine", "Histidine Positive",
          "Cystein SSbond", "Aspartic acid Neutral","Glutatmic acid Neutral", "Lysine Neutral", "Cystein Negative",

          "phosphothreonine", "phosphoserine","o-phosphotyrosine","n-dimethyl-lysine","selenomethionine", "beta-aspartyl",
          "alpha-aminobutyric acid", "s-hydroxycysteine", "s-cysteinesulfinic acid","s,s-(2-hydroxyethyl)thiocysteine", "s-oxy cysteine", "o-methylcysteine",
          "s-oxymethionine", "4-fluoro-l-phenylalanine", "lysine nz-carboxylic acid", "cysteine-s-dioxide", "cysteinesulfonic acid" )

## Formula of aa (from HICUP - e.g. http://xray.bmc.uu.se/hicup/CSW/)
formula <- c("C3 H7 N O2",    "C6 H15 N4 O2",  "C4 H8 N2 O3",   "C4 H6 N O4",     "C3 H7 N O2 S",
             "C5 H10 N2 O3",  "C5 H8 N O4",    "C2 H5 N O2",    "C6 H9 N3 O2",    "C6 H13 N O2",
             "C6 H13 N O2",   "C6 H15 N2 O2",  "C5 H11 N O2 S", "C9 H11 N O2",    "C5 H9 N O2",
             "C3 H7 N O3",    "C4 H9 N O3",    "C11 H12 N2 O2", "C9 H11 N O3",    "C5 H11 N O2",
             
             "C6 H9 N3 O2",   "C6 H9 N3 O2",   "C6 H10 N3 O2",  "C6 H9 N3 O2",    "C6 H9 N3 O2",   "C6 H10 N3 O2",
             "C3 H6 N O2 S",  "C4 H8 N O4",    "C5 H9 N O4",    "C6 H14 N2 O2",   "C3 H6 N O2 S",
             
             "C4 H10 N O6 P", "C3 H8 N O6 P",  "C9 H12 N O6 P", "C8 H18 N2 O2",   "C5 H11 N O2 SE", "C4 H7 N O4",
             "C4 H9 N1 O2",   "C3 H7 N O3 S",  "C3 H7 N O4 S",  "C5 H11 N O3 S2", "C3 H7 N O3 S",   "C4 H9 N O2 S",
             "C5 H11 N O3 S", "C9 H10 F N O2", "C7 H14 N2 O4",  "C3 H7 N O4 S",   "C3 H7 N O5 S")

## Formula of peptide aa residue
formula <- c("C3 H5 N O1",    "C6 H13 N4 O1",  "C4 H6 N2 O2",   "C4 H4 N O3",    "C3 H5 N O1 S",
             "C4 H9 N2 O2",   "C5 H6 N O3",    "C2 H3 N O1",    "C6 H7 N3 O1",   "C6 H11 N O1",
             "C6 H11 N O1",   "C6 H13 N2 O1",  "C5 H9 N O1 S",  "C9 H9 N O1",    "C5 H7 N O1",
             "C3 H5 N O2",    "C4 H7 N O2",    "C11 H10 N2 O1", "C9 H9 N O2",    "C5 H9 N O1",
             
             "C6 H7 N3 O1",   "C6 H7 N3 O1",   "C6 H8 N3 O1",   "C6 H7 N3 O1",    "C6 H7 N3 O1",   "C6 H8 N3 O1",
             "C3 H4 N O1 S",  "C4 H5 N O3",    "C5 H7 N O3",    "C6 H12 N2 O1",   "C3 H4 N O1 S",
             
             "C4 H8 N O5 P",  "C3 H6 N O5 P",  "C9 H10 N O5 P", "C8 H16 N2 O1",   "C5 H9 N O1 SE", "C4 H5 N O3",
             "C4 H7 N1 O1",   "C3 H5 N O2 S",  "C3 H3 N O3 S",  "C5 H9 N O2 S2",  "C3 H5 N O2 S",  "C4 H7 N O1 S",
             "C5 H9 N O2 S",  "C9 H8 F N O1",  "C7 H12 N2 O3",  "C3 H5 N O3 S",   "C3 H7 N O5 S")


## Build the matrix
mat <- cbind(aa, aa321(aa), formula, desc)

## calculate the masses from the formulas
pept <- apply(mat, 1, function(x) formula2mass(x["formula"]))

mat <- cbind(mat, pept)
mat <- mat[,c(1,2,5,3,4)]

## standard residues
mat.st <- mat[1:20,]

## order the non-standard residues
inds <- seq(21, nrow(mat))
mat.spec <- mat[inds,]

inds.o <- order(mat.spec[,"aa"])
mat.spec <- mat.spec[inds.o,]

## column names
new.mat <- rbind(mat.st, mat.spec)
colnames(new.mat) <- c("aa3", "aa1", "aaMass", "name", "formula")
rownames(new.mat) <- new.mat[,"aa3"]

## write matrix to file
write.table(new.mat, quote=T, file="../bio3d/inst/matrices/aa_mass.mat.new", row.names=TRUE)


