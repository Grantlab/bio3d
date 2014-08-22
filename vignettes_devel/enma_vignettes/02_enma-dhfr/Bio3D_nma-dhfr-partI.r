#' # Supporting Material: Integrated structural and evolutionary ensemble analysis with Bio3D
#' ## Lars Skjaerven, Xin-Qiu Yao & Barry J. Grant

#+ setup, include=FALSE
opts_chunk$set(dev='pdf')

#+ preamble, include=FALSE, eval=FALSE
library(knitr)
spin('Bio3D_nma-dhfr-partI.r')
system("pandoc -o Bio3D_nma-dhfr-partI.pdf Bio3D_nma-dhfr-partI.md")

#' ## Background:
#' Bio3D is an R package that provides interactive tools for structural bioinformatics. The primary focus of Bio3D is the analysis of bimolecular structure, sequence and simulation data.

#'
#' #### Requirements:
#' Detailed instructions for obtaining and installing the Bio3D package on various platforms can be found in the [**Installing Bio3D Vignette**](http://thegrantlab.org/bio3d/download/download.html) available both on-line and from within the Bio3D package. In addition to Bio3D the _MUSCLE_ and _CLUSTALO_ multiple sequence alignment programs (available from the [muscle home page](http://www.drive5.com/muscle/) and [clustalo home page](http://www.clustal.org/omega)) must be installed on your system and in the search path for executables. Please see the installation vignette for further details.


#'
#' ## Part I: Ensemble NMA - Bacterial DHFR structures
#' In this vignette we perform *ensemble normal modes analysis* (NMA) on the complete
#' collection of *E.coli* Dihydrofolate reductase (DHFR) structures in the protein databank (PDB).
#' Starting from only one PDB identifier (1rx2) we show how to search using BLAST, fetch and align
#' the structures, and finally calculate the normal modes of each individual structre in order
#' to probe for potential differences in structural flexibility.


#+ example1a, results="hide"
library(bio3d)

#' ### Blast search
#' Below we perform a blast search of the PDB database to identify related structures 
#' to our query *E.coli* DHFR PDB structure. In this particular example we use function 
#' **get.seq(PDBID)** to fetch the query sequence as input to **blast.pdb()**. Note that 
#' **get.seq()** would also allow the UniProt identifier. 

#+ example1b, cache=TRUE, warning=FALSE,
if(file.exists("blast-ecoliDHFR.RData")) {
  load("blast-ecoliDHFR.RData")
} else {
  blast <- blast.pdb(get.seq("1rx2_A"))
  save(blast, file="blast-ecoliDHFR.RData")  
}

#'
#' Function **plot.blast()** facilitates the visualization and filtering of the Blast results. It will attempt to set a seed position to the point of largest drop-off in normalized scores (i.e. the biggest jump in  E-values). In this particular case we specify a cutoff of 225 to include only the relevant *E.coli* structures:

#+ fig1-1, fig.cap="Blast results. Visualize and filter blast results through function **plot.blast()**. Here we proceed with 90 of the top scoring hits."
#+ example1c, cache=TRUE, warning=FALSE,
hits <- plot(blast, cutoff=225)

#'
#' The blast result identified a total of 90 related PDB structures to our query sequence. 
#' The PDB identifiers of this collection are accessible through the attribute `hits$pdb.id`. 
#' Note that adjusting the cutoff argument (to **plot.blast()**) will result in a decrease or
#' increase of hits. 

#' 
#' We can now use function **get.pdb()** and **pdbslit()** to fetch and parse the identified
#' structures. Finally, we use **pdbaln()** to align the PDB structures. 


#+ example1d1, cache=TRUE, warning=FALSE, results='hide'
# fetch PDBs
raw.files <- get.pdb(hits$pdb.id, path = "raw_pdbs")

#+ example1d2, cache=TRUE, warning=FALSE, results='hide'
# split by chain ID
files <- pdbsplit(raw.files, ids = hits$pdb.id, path = "raw_pdbs/split_chain", ncore=4)

#+ example1d3, cache=TRUE, warning=FALSE, results='hide'
# align structures
pdbs.all <- pdbaln(files, fit=TRUE)

#'
#' The *pdbs* object now contains *aligned* C-alpha atom data, including Cartesian coordinates,
#' residue numbers, residue types, and B-factors. The sequence alignment is also stored by default
#' to the FASTA format file 'aln.fa' (to view this you can use an alignment viewer such as SEAVIEW,
#' see _Requirements_ section above). 

#'
#' At this point the *pdbs* object contains all identified 90 structures. This might include 
#' structures with missing residues, and/or multiple structurally redundant conformers. For our
#' subsequent NMA missing in-structure residues might bias the calculation, and redundant structures
#' can be useful to omit to reduce the computational load. Below we inspect the connectivity of the
#' PDB structures with a function call to **inspect.connectivity()**, and **pdbs.filter()** to
#' filter out those structures from our *pdbs* object. Similarly, we  omit structures that are
#' conformationally redundant to reduce the computational load with funtion **rmsd.filter()**:

#+ example1e, cache=TRUE, warning=FALSE, results='hide'
# remove structures with missing residues
conn <- inspect.connectivity(pdbs.all, cut=4.05)
pdbs <- pdbs.filter(pdbs.all, row.inds=which(conn))

# which structures are omitted
which(!conn)

# remove conformational redundant structures
rd <- rmsd.filter(pdbs$xyz, cutoff=0.1, fit=TRUE)
pdbs <- pdbs.filter(pdbs, row.inds=rd$ind)

# remove "humanized" e-coli dhfr
excl <- unlist(lapply(c("3QL0", "4GH8"), grep, pdbs$id))
pdbs <- pdbs.filter(pdbs, row.inds=which(!(1:length(pdbs$id) %in% excl)))

# a list of PDB codes of our final selection
ids <- unlist(strsplit(basename(pdbs$id), split=".pdb"))

#'
#' Use **print()** to see a short summary of the pdbs object:
#+ example1f, warning=FALSE, cache=TRUE
print(pdbs, alignment=FALSE)


#'
#' ### Annotate collected PDB structures 
#' Function **pdb.annotate()** provides a convenient way of annotating the PDB
#' files we have collected. Below we use the function to annotate each structure to its
#' source species. This will come in handy when annotating plots later on:

#+ example1g, warning=FALSE, cache=TRUE, results='hide'
anno <- pdb.annotate(ids)
print(unique(anno$source))


#'
#' ### Principal component analysis
#' A principal component analysis (PCA) can be performed on the structural ensemble
#' (stored in the *pdbs* object) with function **pca.xyz()**. To obtain meaningful results
#' we first superimpose all structures on the *invariant core* (function **core.find()**).

#+ example1h, warning=FALSE, cache=TRUE, results='hide'
core <- core.find(pdbs)
pdbs$xyz = pdbfit(pdbs, core$c0.5A.xyz)

gaps.pos <- gap.inspect(pdbs$xyz)
gaps.res <- gap.inspect(pdbs$ali)
pc.xray <- pca.xyz(pdbs$xyz[,gaps.pos$f.inds])

mktrj(pc.xray, pc=1,
      resno=pdbs$resno[1, gaps.res$f.inds],
      resid=pdbs$resid[1, gaps.res$f.inds])

#'
#' Function **rmsd()** are then used to calculate all pairwise RMSD values in order to
#' perform a quick clustering analysis:

#+ example1i, warning=FALSE, cache=TRUE, results='hide'
rd <- rmsd(pdbs)
hc.rd <- hclust(dist(rd))
grps.rd <- cutree(hc.rd, k=4)


#+ figi-1, fig.cap="Projection of MD conformers onto the X-ray PC space."
plot(pc.xray$z[,1:2], col="grey50", pch=16, cex=1.3,
     ylab="Prinipcal Component 2", xlab="Principal Component 1")
points(pc.xray$z[,1:2], col=grps.rd, pch=16, cex=0.9)


#+ example1j, eval=FALSE
pdbfit(pdbs.filter(pdbs, row.inds=which(grps.rd==1)), outpath="grps1") ## closed
pdbfit(pdbs.filter(pdbs, row.inds=which(grps.rd==2)), outpath="grps2") ## open
pdbfit(pdbs.filter(pdbs, row.inds=which(grps.rd==3)), outpath="grps3") ## occluded


#'
#' ### Normal modes analysis
#' Function **nma.pdbs()** will calculate the normal modes of each protein structures stored
#' in the *pdbs* object. The normal modes are calculated on the full structures as provided
#' by object *pdbs*. With the default argument `rm.gaps=TRUE` unaligned atoms 
#' are omitted from output:

#+ example1k, cache=TRUE, warning=FALSE,
modes <- nma.pdbs(pdbs, rm.gaps=TRUE, ncore=4)

#'
#' The *modes* object of class *enma* contains aligned normal mode data including fluctuations,
#' RMSIP data, and aligned eigenvectors. A short summary of the *modes* object can be obtain by
#' calling the function **print()**, and the aligned fluctuations can be plotted with function
#' **plot.enma()**. 

#+ example1l, cache=TRUE, warning=FALSE,
print(modes)

#+ fig1l-1, fig.cap="Normal mode fluctuations of the E.coli ensemble", fig.height=4.5,
# plot modes fluctuations
plot(modes, pdbs=pdbs, col=grps.rd)

#' **plot.enma()** also facilitates the calculation and visualization of sequence conservation
#' (use argument `conservation=TRUE`), and fluctuation variance (use argument `variance=TRUE`).

#+ example1m, cache=TRUE, warning=FALSE,
hc.nma <- hclust(as.dist(1-modes$rmsip))
grps.nma <- cutree(hc.nma, k=4)

#+ fig1m-1, fig.cap="Clustering of structures based on their normal modes similarity"
heatmap(1-modes$rmsip, distfun = as.dist, labRow = ids, labCol = ids,
        ColSideColors=as.character(grps.nma), RowSideColors=as.character(grps.nma))

#'
#' ### Fluctuation analysis
#' Comparing the mode fluctuations of two groups of structures 
#' can reveal specific regions of distinc flexibility patterns. Below we focus on the
#' differences between the open (black), closed (red) and occluded (green) conformations
#' of the *E.coli* structures: 


           
#+ fig1n-1, fig.cap="Comparison of mode fluctuations between open (black) and closed (red).", fig.height=4.5,
cols <- grps.rd
cols[which(cols!=1 & cols!=2)]=NA
plot(modes, pdbs=pdbs, col=cols, signif=TRUE)

#+ fig1n-2, fig.cap="... open (black) and occluded (green).", fig.height=4.5,
cols <- grps.rd
cols[which(cols!=1 & cols!=3)]=NA
plot(modes, pdbs=pdbs, col=cols, signif=TRUE)

#+ fig1n-3, fig.cap="... closed (red) and occluded (green).", fig.height=4.5,
cols <- grps.rd
cols[which(grps.rd!=2 & grps.rd!=3)]=NA
plot(modes, pdbs=pdbs, col=cols, signif=TRUE)





#'
#' ### Compare with MD simulation
#+ example1o, eval=TRUE, results='hide'
mddir <- "../01_md-simulations/05_1rx2"
pdbmd <- paste(mddir, "01_build", "prot.pdb", sep="/")

pdb <- read.pdb(pdbmd)
ca.inds <- atom.select(pdb, 'calpha')
pdb <- trim.pdb(pdb, ca.inds)

trj <- read.ncdf(paste(mddir, "03_prod", "prod_noWAT.nc", sep="/"), at.sel=ca.inds)
md.inds <- pdb2aln.ind(aln=pdbs, pdb=pdb, id="md", inds=gaps.res$f.inds)

trj=trj[, atom2xyz(md.inds)]
trj=fit.xyz(fixed=pdbs$xyz[1, ],  mobile=trj,
    fixed.inds=core$c0.5A.xyz, mobile.inds=core$c0.5A.xyz, outpath = "fitlsq-md/")

proj <- pca.project(trj, pc.xray)
cols <- densCols(proj[,1:2])

#+ fig1o-1, fig.cap="Projection of MD conformers onto the X-ray PC space."
plot(proj[,1:2], col=cols, pch=16,
     ylab="Prinipcal Component 2", xlab="Principal Component 1",
     xlim=range(pc.xray$z[,1]), ylim=range(pc.xray$z[,2]))
points(pc.xray$z[,1:2], col=1, pch=1, cex=1.1)
points(pc.xray$z[,1:2], col=grps.rd, pch=16)

#+ example1p, eval=TRUE
# PCA of the MD trajectory
pc.md <- pca.xyz(trj)

# compare MD-PCA and NMA
r <- rmsip(pc.md$U, modes$U.subspace[,,1])

print(r)

#+ fig1p-1, fig.cap="RMSIP map betwen normal modes and principal components of a 5 ns long MD simulation."
plot(r, xlab="MD PCA", ylab="NMA")


#' ## Document Details
#' This document is shipped with the Bio3D package in both R and PDF formats. All code can be extracted and automatically executed to generate Figures and/or the PDF with the following commands:

#+ close, include=TRUE, eval=FALSE
library(knitr)
spin('Bio3D_nma-dhfr-part0.r')
system("pandoc -o Bio3D_nma-dhfr-partI.pdf Bio3D_nma-dhfr-partI.md")

#'
#' ## Information About the Current Bio3D Session
#'
print(sessionInfo(), FALSE)





