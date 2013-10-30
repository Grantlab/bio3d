#' # Supporting Material: Enhanced Methods for Normal Mode Analysis with Bio3D
#' ## Lars Skjaerven, Xin-Qiu Yao & Barry J. Grant

#+ preamble, include=FALSE, eval=FALSE
library(knitr)
spin('Bio3D_nma.r')
system("pandoc -o Bio3D_nma.pdf Bio3D_nma.md")


#' ## Background:
#' Bio3D[^1] is an R package that provides interactive tools for structural bioinformatics. The primary focus of Bio3D is the analysis of bimolecular structure, sequence and simulation data.
#'
#' Normal mode analysis (NMA) is one of the major simulation techniques used to probe large-scale motions in biomolecules. Typical application is for the prediction of functional motions in proteins. Version 2.0 of the Bio3D package now includes extensive NMA facilities. These include a unique collection of multiple elastic network model force-fields (see Example 1 below), automated ensemble analysis methods (Example 2) and variance weighted NMA (Example 3). Here we demonstrate the use of these new features with working code that comprise complete executable examples[^2].
#'
#' [^1]: The latest version of the package, full documentation and further vignettes (including detailed installation instructions) can be obtained from the main Bio3D website: http://thegrantlab.org/bio3d/
#'
#' [^2]: This document contains executable code that generates all figures contained within this document. See help(vignette) within R for further details. 
#'

#'
#' ## Example 1A: Basic normal mode analysis
#' Normal mode analysis (NMA) of a single protein structure can be carried out by providing a PDB object to the function **nma()**. In the code below we first load the Bio3D package and then download an example structure of hen egg white lysozyme (PDB id *1hel*) with the function **read.pdb()**. Finally the function **nma()** is used perform the normal mode calculation:

#+ example1_A, results="hide"
library(bio3d)
pdb <- read.pdb("1hel")
modes <- nma(pdb)

#' A short summary of the returned *nma* object contained within the new variable *modes* can be obtained by simply calling the function **print()**:

print(modes)

#' This reveals the function call resulting in the *nma* object along with the total number of stored normal modes. For PDB id *1hel* there are 129 amino acid residues, and thus $3*129=387$ modes in this object (in which the first six are so-called trivial modes with zero frequency corresponding to rigid-body rotation and translation). The frequency of the next six lowest-frequency modes is also printed. Note that the returned *nma* object consists of a number of attributes listed on the _+attr:_ line. These attributes contain the detailed results of the calculation and a complete description of their contents and structure can be found on the **nma()** functions help page accessible with the command: `help(nma)`. Additional examples of data access are provided further below. However, to get a quick overview of the results one can simply call the **plot()** function on the returned *nma* object. This will produce a summary plot of (1) the eigenvalues, (2) the mode frequencies, and (3) the atomic fluctuations (See Figure 1).

#+ fig1_a_1, fig.cap="Summary plot of NMA results for hen egg white lysozyme (PDB id *1hel*). The optional *sse=pdb* argument provided to **plot.nma()** results in a secondary structure schematic being added to the top and bottom margins of the fluctuation plot (helices black and strands gray). Note the larger fluctuations predicted for loop regions." 
plot(modes, sse=pdb)


#'
#' ## Example 1B: Multiple force fields for normal mode analysis
#' The main Bio3D normal mode analysis function, **nma()**, requires a set of coordinates, as obtained from the **read.pdb()** function, and the specification of a force field describing the interactions between constituent atoms. By default the *calpha* force field originally developed by Konrad Hinsen is utilized. This employs a spring force constant differentiating between nearest-neighbor pairs along the backbone and all other pairs. The force constant function was parameterized by fitting to a local minimum of a crambin model using the AMBER94 force field. However, a large number of additional force fields are also available. Full details of these force fields can be obtained with the command `help(load.enmff)`. With the code below we briefly demonstrate their usage and comparison: 


#+ help, eval=FALSE
help(load.enmff)

#+ example1_B, results="hide"
modes.a <- nma(pdb, ff="calpha")
modes.b <- nma(pdb, ff="anm")
modes.c <- nma(pdb, ff="pfanm")
modes.d <- nma(pdb, ff="calphax")
modes.e <- nma(pdb, ff="reach")

#+ example1_B_rmsip
rmsip(modes.a, modes.b)
## RMSIP matrix plot here with adjplot or similar...
#adjplot( rmsip(modes.a, modes.b, modes.c, modes.d, modes.e) )

#' Compare to expermental difference vector here also and generate molecular figure

#'
#' ## Example 2: Ensemble normal mode analysis
#' The analysis of multiple protein structures (e.g. a protein family) can be accomplished with the **nma.pdbs()** function. This will take aligned input structures, as generated by the **pdbaln()** function for example, and perform NMA on each structure collecting the results in manner that facilitates the interoperation of similarity and dissimilarity trends in the structure set. Here we will analyze a collection of protein kinase structures with low sequence identity (2A) and large set of closely related transducin heterotrimeric G protein family members (2B). 

#'
#' ### Example 2A: Protein kinases 

#+ example2_A, cache=TRUE, results="hide"
# Select Protein Kinase PDB IDs
ids <- c("4b7t_A", "2exm_A", "1opj_A", 
         "4jaj_A", "1a9u_A", "1tki_A", 
          "1phk_A", "1csn_A", "1lp4_A") 
         # "1m14", "3dnd", "2jdo", "1JKL", "1gng", "2src", "1OMW", "1b6c_D",

raw.files <- get.pdb(ids, path="raw_pdbs")
files <- pdbsplit( raw.files, ids )

# Alignment of structures
pdbs <- pdbaln(files)

# NMA on all structures
modes <- nma.pdbs(pdbs, full=TRUE)

#+ plot_enma, width=7, height=4, fig.cap="Results of ensemble NMA on select protein kinase superfamily members"
# Plot fluctuation data
plot(modes, typ="l")
## add legend with PDB IDs and cols


#'
#' ### Example 2B: Transducin
#' 
#' This example will run **nma.pdbs()** on the transducin family containing 53 PDB structures.

#' First, we download pre-compiled structure and alignment data for transducin.
#+ data, cache=TRUE
# We may think of puting this data (or a subset) within the package!
download.file("http://www-personal.umich.edu/~xinqyao/transducin.RData", 
              "transducin.RData")
load("transducin.RData")

#' The object **pdbs** contains the C-alpha atoms coordinates of the 53 structures with
#' atomic positions aligned based on the multiple sequence alignment. 
#+, cache=TRUE, results="hide"
nmodes <- nma.pdbs(pdbs)
ind <- grep("1TAG_A", pdbs$id)
pdb <- read.pdb("1tag")
resno <- pdbs$resno[ind, ]
resno[!is.na(resno)] <- pdb$atom[pdb$calpha,"resno"]
sse<-dssp(pdb)

#' Objects **gaps.atom** and **gaps.xyz** are prepared previously with calling
#' gap.inspect(pdbs\$ali) and gap.inspect(pdbs\$xyz), respectively. **ligs** annotates
#' the nucleotide state of each structure based on a preliminary inspection on the structure data, 
#' while **vcolors** accordingly colors nucleotide state with green (GDP) and red (GTP).
#'

#+ example2_B, fig.cap="Atomic fluctuation of transducin predicted by NMA"
plot(nmodes, col=vcolors)
 
##  plot(resno[gaps.atom$f.inds], nmodes$fluctuations[i,], type='h', col=vcolors[i], sse=sse)

text(x=176, y=1, label="SW I")
text(x=200, y=1.3, label="SW II")
text(x=210, y=2.2, label="SW III")


#' The inter-structure relationship can be characterized via comparing the modes predicted by NMA. 
#' The similarity of structures in terms of dynamic property is calculated with the root mean square
#' inner product (RMSIP) of low-frequency vibrational modes. As a comparison, we also calculated the
#' inter-structure root mean square deviation (RMSD) of the C-alpha atoms.
#+ example2_B2, fig.cap="Distribution of overlap value among transducin family (Do we need this figure?)"
rmsip.map <- nmodes$rmsip
rownames(rmsip.map) <- substr(basename(pdbs$id), 1, 6)
colnames(rmsip.map) <- substr(basename(pdbs$id), 1, 6)
hist(rmsip.map, breaks=70)
#+ example2_B3, fig.cap="RMSIP matrix of transducin family"
heatmap((1-rmsip.map), labRow=ligs, symm=TRUE)

#+, cache=TRUE
rmsd.map <- rmsd(pdbs$xyz, a.inds=gaps.xyz$f.inds, fit=TRUE)
dimnames(rmsd.map) <- dimnames(rmsip.map)
#+ example2_B4, fig.cap="RMSD matrix of transducin family"
heatmap(rmsd.map, labRow=ligs, symm=TRUE)

#'
#' ## Example 3: Variance weighted normal mode analysis
#' In this example we illustrate an approach of weighting the pair force constants based on the variance of the inter atomic distances obtained from an ensemble of structures (e.g. available X-ray structures). The motivation for such variance-weighting is to reduce the well known depence of the force constants on the one structure upon which they derived derived [ADD REFS here!].
#'
#' ### Example 3A: GroEL
#' The GroEL subunit consists of 524 residues which comprise three distinct domains inter-connected by two hinge regions facilitating large conformational changes. In the example below we will see how the normal modes of the open state corresponds nicely with the observed conformational change (by X-ray and EM studies), contrary to the closed state. We will then use an ensemble of X-ray structures as weighting to the pair-force constants:

#+, cache=TRUE, results="hide"
# Set the ensemble PDB ids
ids <- c("1sx4_[A-C]", "1xck_[A-C]")

##, "1svt_[A-C]", "1sx3_[A-C]", "1kp8_[A-C]")

# Download and split PDBs by chain ID
raw.files <- get.pdb(ids, path = "raw_pdbs", gzip=TRUE)
files <- pdbsplit(raw.files, ids, path = "raw_pdbs/split_chain/")

# Align and superimpose coordinates
pdbs <- pdbaln(files, fit=TRUE)

#'
#' The 'pdbs' object contains *aligned* C-alpha atom data, including cartesian coordinates, residue numbers, residue type and B-factors. The sequence alignement is also stored by default to the fasta file 'aln.fa'.
#'

#'
#' #### Calculate normal modes
#' Next we will calculate the normal modes of the open and closed conformational state. They are stored at indices 1 and 4, respectively, in our 'pdbs' object. Use the **pdbs2pdb()** to fetch the pdb objects which is needed for the input to **nma()**. 
#'

#+, cache=TRUE,
## Inspect gaps
gaps.res <- gap.inspect(pdbs$ali)
gaps.pos <- gap.inspect(pdbs$xyz)

## Convert back to access PDB objects
pdb.list   <- pdbs2pdb(pdbs, inds=c(1,4))

pdb.open   <- pdb.list[["1sx4_A"]]
pdb.closed <- pdb.list[["1xck_A"]]

modes.open   <- nma(pdb.open)
modes.closed <- nma(pdb.closed)

#'
#' #### Overlap analysis
#' Use overlap analysis to determine the agreement between the normal mode vectors and the conformational difference vector:
#'

#+, cache=TRUE,
## Difference vector
diff.vec <- difference.vector(pdbs$xyz[c(1,4), gaps.pos$f.inds])

## Calculate overlap
oa <- overlap(modes.open,   diff.vec)
ob <- overlap(modes.closed, diff.vec)


#+ example3A_overlapplotA, fig.width=7, fig.height=6, fig.cap="Overlap anlaysis"
plot(oa$overlap.cum[1:10], type='b', ylim=c(0,1),
     ylab="Squared overlap", xlab="Mode index",
     cex.lab=1.4, axes=F, lwd=2)
lines(ob$overlap.cum[1:10], type='b', lty=2, col=2, lwd=2)

legend("bottomright",
       c("Open state", "Closed state"),
       col=c("black", "red"), lty=c(1,2))

box()
axis(1, cex.axis=1.3)
axis(2, cex.axis=1.3)

#'
#' #### Variance weighting
#' From the overlap analysis above we see the remarkable agreement between the conformational difference vector and the normal modes calculated on the open structures. Contrary, the lowest frequency modes of the closed structures does not show the same behaviour. We will thus proceed with the weighting of the force constants. First we'll define a quick function for calculating the weights which takes a matrix of cartesian coordinates as input:
#'

#+, cache=TRUE,
"make.weights" <- function(xyz) {
  ## Calculate pairwise distances
  natoms <- ncol(xyz) / 3
  all <- array(0, dim=c(natoms,natoms,nrow(xyz)))
  for( i in 1:nrow(xyz) ) {
    dists <- dist.xyz(xyz[i,])
    all[,,i] <- dists
  }
  
  ## Calculate variance of pairwise distances
  all.vars <- apply(all, 1:2, var)
  
  ## Make the final weights
  weights <- 1 - (all.vars / max(all.vars))
  return(weights)
}

#+, cache=TRUE
## Calcualte the weights
wts <- make.weights(pdbs$xyz[, gaps.pos$f.inds])

#' Weights to the force constants can be included by the argument 'fc.weights' to function **nma()**. 
#' This needs be a matrix with dimensions NxN. Here we will run a small for-loop with increasing the 
#' strength of the weigthing at each step and store the new overlap values in the variable 'ocs':

#+, cache=TRUE, results="hide"
ocs <- NULL
for ( i in 1:10 ) {
  modes.wtd <- nma(pdb.closed, fc.weights=wts**i)
  oc        <- overlap(modes.wtd,    diff.vec)
  ocs       <- rbind(ocs, oc$overlap.cum)
}

#+ example3A_overlapplot, fig.width=7, fig.height=6, fig.cap="Overlap plot"
plot(oa$overlap.cum[1:10], type='b', ylim=c(0,1),
     ylab="Squared overlap", xlab="Mode index",
     cex.lab=1.4, axes=F, lwd=2)
lines(ob$overlap.cum[1:10], type='b', lty=2, col=1, lwd=2)

cols <- rainbow(10)
for ( i in 1:nrow(ocs) ) {
  lines(ocs[i,1:10], type='b', lty=1, col=cols[i])
}

legend("bottomright",
       c("Open state", "Closed state", "Closed state (weighted)"),
       col=c("black", "black", "green"), lty=c(1,2,1))

box()
axis(1, cex.axis=1.3)
axis(2, cex.axis=1.3)


#'
#' #### RMSIP calculation
#' Root mean square inner product (RMSIP) can be used to compare mode subspaces:
#+ example3A_rmsip
ra <- rmsip(modes.open, modes.wtd)
rb <- rmsip(modes.open, modes.closed)

ra$rmsip
rb$rmsip

#+ example3A_rmsipplot, fig.width=9, fig.height=4, fig.cap="RMSIP map"
par(mfrow=c(1,2))
image(1:10, 1:10, ra$overlap, col=gray(50:0/50), zlim=c(0,1),
      ylab="NMA(open)", xlab="NMA(weighted)")
image(1:10, 1:10, rb$overlap, col=gray(50:0/50), zlim=c(0,1),
      ylab="NMA(open)", xlab="NMA(closed)")


#'
#' #### Match with PCA
#'

## Calculate the PCs
pc.xray <- pca.xyz(pdbs$xyz[,gaps.pos$f.inds])

rmsip(pc.xray, modes.closed)
rmsip(pc.xray, modes.wtd)








#'
#' ### Example 3B: Transducin
#' 
#' This example will run **nma()** on transducin with various variance weighted force constants. 
#' The modes predicted by NMA will be compared with PCA results over the transducin family.

#' First, download data and wrap a function for making variance weights:
#+, cache=TRUE
# We may think of puting the online data somewhere else!!!
download.file("http://www-personal.umich.edu/~xinqyao/transducin.RData", 
              "transducin.RData")
load("transducin.RData")

mk.weights <- function(xyz) {
  # weights <- mk.weights(xyz)
  natoms <- ncol(xyz) / 3
  all <- array(0, dim=c(natoms,natoms,nrow(xyz)))
  for( i in 1:nrow(xyz) ) {
    dists <- dist.xyz(xyz[i,])
    all[,,i] <- dists
  }
  all.vars <- apply(all, 1:2, var)
  cat( paste("Max:", round(max(all.vars),2),
             "  Min:", round(min(all.vars),2) ),"\n" )
  return(1 - (all.vars / max(all.vars)))
}

#' Then, we take one structure for each GDP and GTP state, and run **nma()** with force field **calpha**.
#' The PDB objects **gdp** and **gtp** contains the structural information for the C-alpha atoms of 
#' transducin in GDP (PDB ID 1TAG) and GTP (PDB ID 1TND) states, respectively. The coordinates 
#' in **gdp** and **gtp** were fitted to the **pdbs** object based on all non-gap C-alpha positions. 
#+, results="hide"
inds <- sapply(c("1TAG_A", "1TND_B"), grep, pdbs$id)
inds.gdp <- atom.select(gdp, resno=pdbs$resno[inds[1], gaps.atom$f.inds])
inds.gtp <- atom.select(gtp, resno=pdbs$resno[inds[2], gaps.atom$f.inds])

#+
modes.gdp <- nma(gdp, inds=inds.gdp, fc.weights=NULL)
modes.gtp <- nma(gtp, inds=inds.gtp, fc.weights=NULL)

#' Now, we calculate the pairwise distance variance based on the structure ensemble 
#' with the wrapped function mentioned above. This will be used to weight the force constants
#' in the elastic network model. The object **xyz** is a numeric matrix containing 
#' aligned cartesian coordinates all fitted to the first structure in **pdbs**, 
#' xyz <- pdbfit(pdbs).
weights <- mk.weights(xyz[, gaps.xyz$f.inds])

modes.gdp.b <- nma(gdp, inds=inds.gdp, fc.weights=weights**100)
modes.gtp.b <- nma(gtp, inds=inds.gtp, fc.weights=weights**100)

#' To evaluate the results, we calculate the overlap (square dot product) between modes predicted by 
#' variance weighted or non-weighted NMA and the first principle component from PCA (contained
#' in the object **pc.xray** returned by calling **pca.xyz(xyz[, gaps.xyz$f.inds])**.
oa <- overlap(modes.gdp, pc.xray$U[,1])
ob <- overlap(modes.gtp, pc.xray$U[,1])
oc <- overlap(modes.gdp.b, pc.xray$U[,1])
od <- overlap(modes.gtp.b, pc.xray$U[,1])

#+, fig.cap="Variance weighted force constants improve NMA prediction"
plot(oa$overlap.cum, type='o', ylim=c(0,1), col="darkgreen", lwd=2, xlab="Mode", 
     ylab="Cummulative overlap")
lines(ob$overlap.cum, type='o', ylim=c(0,1), col="red", lwd=2)
lines(oc$overlap.cum, type='b', ylim=c(0,1), col="darkgreen", lwd=2)
lines(od$overlap.cum, type='b', ylim=c(0,1), col="red", lwd=2)
text(20, oa$overlap.cum[20], label=round(oa$overlap.cum[20], 2), pos=3)
text(20, ob$overlap.cum[20], label=round(ob$overlap.cum[20], 2), pos=3)
text(20, oc$overlap.cum[20], label=round(oc$overlap.cum[20], 2), pos=3)
text(20, od$overlap.cum[20], label=round(od$overlap.cum[20], 2), pos=3)
legend(x=2, y=1.0, pch=1, lty=c(1, 1, 2, 2), col=c("darkgreen", "red", 
       "darkgreen", "red"), legend=c("GDP", "GTP", "Weighted GDP", "Weighted GTP"))

#'
#' ## Document Details
#' This document is shipped with the Bio3D package in both latex and PDF formats. All code can be extracted and automatically executed to generate Figures and/or PDF with the following commands:

#+ close, include=TRUE, eval=FALSE
library(knitr)
spin('Bio3D_nma.r')
system("pandoc -o Bio3D_nma.pdf Bio3D_nma.md")

#'
#' ## Information About the Current Bio3D Session
#'
print(sessionInfo(), FALSE)
