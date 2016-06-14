
## Generate distance restraints for Amber MD run

## Usage:
# prmtop <- read.prmtop("prot.prmtop")
# crd <- read.crd("prot.inpcrd")
# pdb <- as.pdb(prmtop, crd)

## or read pdb generated from
# system("ambpdb -p prot.prmtop < prot.inpcrd > prot.amb.pdb")
# pdb <- read.pdb("prot.amb.pdb")

## Inter-domain restraints:
# d1 <- atom.select(pdb, "noh", resno=8:94)
# d2 <- atom.select(pdb, "noh", resno=1301:1387)
# rst <- amber.rst(pdb, a.inds=d1, b.inds=d2)
# write.rst(rst, file="R.rst")

## Intra-domain restraints:
# d1 <- atom.select(pdb, "noh", resno=8:94)
# rst <- amber.rst(pdb, a.inds=d1)
# write.rst(rst, file="R.rst")

## cut: lower and upper cutoff. look for distances in range (lower, upper]. 
## r: corresponds to r1, r2, r3, r4
## rk: corresponds to rk2, and rk3


amber.rst <- function(pdb, a.inds=NULL, b.inds=a.inds,
                      cut=c(4, 5), r=c(-1.0, -0.25, 0.25, 1.0),
                      rk=c(5, 5), file="dist.rst") {
  
  if(!is.pdb(pdb))
    stop("must supply a pdb or prmtop object")
  
  if(!is.select(a.inds) | !is.select(b.inds))
    stop("must supply atom.selection")

  if(length(cut)!=2)
    stop("cut must be a numeric vector of length 2")

  if(length(r)!=4)
    stop("r must be a numeric vector of length 4")

  if(length(rk)!=2)
    stop("rk must be a numeric vector of length 2")

  if(cut[1]>cut[2])
    stop("cut: lower cutoff must be smaller than upper cutoff")

  if(any(duplicated(pdb$atom$eleno))) {
    stop(paste("Duplicate element numbers (pdb$atom$eleno) found. Re-number, and try again",
               "   e.g. pdb$atom$eleno = 1:length(pdb$atom$eleno) ", sep="\n"))
  }
  
  r <- round(r, 2)
  rk <- round(rk, 2)

  ## distance matrix
  dm <- dist.xyz(pdb$xyz[1, a.inds$xyz], pdb$xyz[1, b.inds$xyz])

  ## matrix names
  res.a <- paste(pdb$atom$resno[a.inds$atom],
                 pdb$atom$chain[a.inds$atom],
                 pdb$atom$insert[a.inds$atom], sep="-")
  res.b <- paste(pdb$atom$resno[b.inds$atom],
                 pdb$atom$chain[b.inds$atom],
                 pdb$atom$insert[b.inds$atom], sep="-")
  
  rownames(dm) <- res.a
  colnames(dm) <- res.b
  #print(dm)
  
  ## avoid duplicate pairs
  if(identical(a.inds, b.inds, ignore.environment = TRUE)) {
    dm[lower.tri(dm, diag=TRUE)] <- NA
  }
  #print(dm)
  
  ## omit covalent atoms
  inds <- which(dm < cut[1])
  dm[inds] <- NA
  #print(dm)


  ## avoid intra-residue restraints
  overlapping <- intersect(res.a, res.b)
  
  if(length(overlapping)>0) {
    for(i in 1:length(unique(overlapping))) {
      tmp.a <- which(res.a==unique(overlapping)[i])
      tmp.b <- which(res.b==unique(overlapping)[i])
      dm[tmp.a, tmp.b] <- NA
    }
  }

  ## finally, collect distances below cutoff
  inds <- which(dm <= cut[2], arr.ind=TRUE)

  if(!length(inds)>0) {
    cat(" no atom distances within provided cutoff range\n")
    return(NULL)
  }

  ## atom numbers
  atoms.a <- as.numeric(pdb$atom$eleno[a.inds$atom][inds[,1]])
  atoms.b <- as.numeric(pdb$atom$eleno[b.inds$atom][inds[,2]])

    
  inds <- which(dm <= cut[2])
  out <- list(atom.a = atoms.a,
              resno.a = pdb$atom$resno[atoms.a],
              elety.a = pdb$atom$elety[atoms.a],
              chain.a = pdb$atom$chain[atoms.a],
              
              atom.b = atoms.b,
              resno.b = pdb$atom$resno[atoms.b],
              elety.b = pdb$atom$elety[atoms.a],
              chain.b = pdb$atom$chain[atoms.b],
              d = round(dm[inds], 2),
              
              r1 = round(dm[inds], 2) + r[1],
              r2 = round(dm[inds], 2) + r[2],
              r3 = round(dm[inds], 2) + r[3],
              r4 = round(dm[inds], 2) + r[4],
              
              rk2 = rk[1],
              rk3 = rk[2])
  
    return(as.data.frame(out, stringsAsFactors=FALSE))
}


write.rst <- function(rst, file="R.rst", append=FALSE) {

    if(!nrow(rst)>0)
        stop("'rst' must contain 1 or more restraints")
    
    rst.lines <- NULL
    for(i in 1:nrow(rst)) {
        
        tmp.lines <- c("#  ",
                       paste("#  Atoms",
                             " ", rst$resno.a[i], "@", rst$elety.a[i],
                             " - ", rst$resno.b[i], "@", rst$elety.b[i],
                             " (", rst$d[i], " Å)",
                             sep=""), 
                       paste(" &rst iat=", rst$atom.a[i], ",", rst$atom.b[i], ", ", sep=""),
                       paste("   r1=", rst$r1[i], ", r2=", rst$r2[i], ", r3=", rst$r3[i], ", r4=", rst$r4[i], ", ", sep=""),
                       paste("   rk2=", rst$rk2[i], ", rk3=", rst$rk3[i], ", &end  ", sep=""))
        
        rst.lines <- c(rst.lines, tmp.lines)
    }
    
    #writeLines(rst.lines, file)
    write(rst.lines, file, append=append)
    cat("  written", nrow(rst), "distance restraints to file:", file, "\n")
    

}
