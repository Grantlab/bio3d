#' # Supplementary Data: Enhanced Methods for Normal Mode Analysis with Bio3D
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
#' In this vignette we calculate normal modes on an ensemble of *E.coli* Dihydrofolate
#' reductase (DHFR) structures, and compare the resulting flexibility profiles with 3 
#' homologous proteins identified through a blast search. We show how structural flexibility
#' relate to sequence conservation, and we identify specific reigions showing
#' distinct fluctuation profiles. 


#+ example1a, results="hide"
library(devtools)
load_all("~/workspace/bio3d/ver_devel/bio3d")

#' ### Blast search
#' Below we perform a blast search of the PDB database to identify related structures 
#' to our query E.coli DHFR PDB structure. In this particular example we use function 
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
#' Visualize and filter the blast results through function **plot.blast()**. The function will attempt to set a seed position to the point of largest drop-off in normalized scores (i.e. the biggest jump in  E-values). In this particular case we specify a cutoff of 92 to include more distant homologues:

#+ example1b1a, eval=FALSE
# default, 91 hits
hits <- plot(blast)

#+ example1b1b, cache=TRUE, warning=FALSE

#+ fig1-0, fig.cap="Blast results. Visualize and filter blast results through function **plot.blast()**. The blast result identified a total of 186 related PDB structures to our query sequence."
# include more hits
hits <- plot(blast, cutoff=92)

#'
#' The blast result identified a total of 186 related PDB structures to our query sequence. 
#' The PDB identifiers of this collection are accessible through the attribute `hits$pdb.id`. 
#' Note that adjusting the 
#' cutoff argument (to **plot.blast()**) will result in a decrease or increase of hits. 

#' 
#' We can now use function **get.pdb()** and **pdbslit()** to fetch and parse the identified
#' structures. Finally, we use **pdbaln()** to align the PDB structures. 

#+ example1b2, cache=TRUE, warning=FALSE, results='hide'
# fetch PDBs
raw.files <- get.pdb(hits$pdb.id, path = "raw_pdbs")

# split by chain ID
files <- pdbsplit(raw.files, ids = hits$pdb.id, path = "raw_pdbs/split_chain", ncore=4)

# align structures
pdbs <- pdbaln(files, fit=TRUE)

#'
#' The *pdbs* object now contains *aligned* C-alpha atom data, including Cartesian coordinates,
#' residue numbers, residue types, and B-factors. The sequence alignment is also stored by default
#' to the FASTA format file 'aln.fa' (to view this you can use an alignment viewer such as SEAVIEW,
#' see _Requirements_ section above). 

#'
#' At this point the *pdbs* object contains all identified 186 structures. This might include 
#' structures with missing residues, and/or multiple structurally redundant conformers. For our
#' subsequent NMA missing in-structure residues might bias the calculation, and redundant structures
#' can be useful to omit to reduce the computational load. Below we inspect the connectivity of the PDB
#' structures with a function call to **inspect.connectivity()**, and **pdbs.filter()** to filter out those
#' structures from our *pdbs* object. Similarly, we  omit structures that are
#' conformationally redundant to reduce the computational load with funtion **rmsd.filter()**:

#+ example1b3, cache=TRUE, warning=FALSE, results='hide'
# remove structures with missing residues
conn <- inspect.connectivity(pdbs, cut=4.05)
pdbs <- pdbs.filter(pdbs, row.inds=which(conn))

# which structures are omitted
which(!conn)

# remove conformational redundant structures
rd <- rmsd.filter(pdbs$xyz, cutoff=0.25, fit=TRUE)
pdbs <- pdbs.filter(pdbs, row.inds=rd$ind)

# a list of PDB codes of our final selection
ids <- unlist(strsplit(basename(pdbs$id), split=".pdb"))

#'
#' Use **print()** to see a short summary of the pdbs object:
#+ example1b4, warning=FALSE, cache=TRUE, eval=TRUE,
print(pdbs, alignment=FALSE)


#'
#' ### Annotate collected PDB structures 
#' Function **pdb.annotate()** provides a convenient way of annotating the PDB
#' files we have collected. Below we use the function to annotate each structure to its
#' source species. This will come in handy when annotating plots later on:

#+ example1c, warning=FALSE, cache=TRUE, 
anno <- pdb.annotate(ids)
print(unique(anno$source))

# Categorize the results
grps <- rep(NA, length(ids))
grps[grep("coli", anno$source)]=1
grps[grep("profunda", anno$source)]=2
grps[grep("anthracis", anno$source)]=3
grps[grep("pestis", anno$source)]=4

mynames <- rep(NA, length(ids))
mynames[grps==1]="E. coli"
mynames[grps==2]="M. profunda"
mynames[grps==3]="B. anthracis"
mynames[grps==4]="Y. pestis"

#'
#' ### Sequence conservation analysis
#' To perform a more rigorous conservation analysis we perform a new blast search,
#' this time using the swissprot database. We align the sequences to our alignment
#' of the PDBs by providing the *pdbs* object as a profile to function **seqaln()**:

#+ example1d, warning=FALSE, cache=TRUE, 
if(file.exists("blast-ecoliDHFR-sp.RData")) {
  load("blast-ecoliDHFR-sp.RData")
} else {
  blast.seqs <- blast.pdb(get.seq("1rx2_A"), database="swissprot")
  save(blast.seqs, file="blast-ecoliDHFR-sp.RData")  
}


#+ fig1-2b, fig.cap="Blast results from the SwissProt database.", fig.height=5,
hits.seq <- plot(blast.seqs, cutoff=60)

#+ example1d2, warning=FALSE, cache=TRUE, results='hide',
seqs <- get.seq(hits.seq$gi.id)
aln <- seqaln(seqs, profile=pdbs, exefile="clustalo", extra.args="--dealign")


#'
#' The conservation for the sequences can be calculated with function **conserv()**.
#' Note that we use only the 27 first entries for the conservation alignment, which is 
#' the number of non-redundant sequences identified from swissprot. 

# there are 27 sequences from swissprot
nseqs <- length(seqs$id)
cons <- conserv(aln$ali[1:nseqs,])

# conservation only for residues in pdbs
gaps.aln <- gap.inspect(aln$ali[(nseqs+1):nrow(aln$ali),])
#cons[gaps.aln$f.inds]

#' 
#' Use function **seqidentity()** to calculate sequence identity for the PDBs:
#+ example1f, warning=FALSE, cache=TRUE, 
seqide <- seqidentity(pdbs)
summary(c(seqide))
hc <- hclust(as.dist(1-seqide))
grps.seq <- cutree(hc, k=4)

#+ fig1-3, fig.cap="Sequence-based clustering of the collected structures shows that E.coli and B.anthracis shows the largest difference in terms of sequence.", fig.height=3.5,
plot(hc, hang=-1, labels=mynames, cex=0.5)


#'
#' ### Normal modes analysis
#' Function **nma.pdbs()** will calculate the normal modes of each protein structures stored
#' in the *pdbs* object. The normal modes are calculated on the full structures as provided
#' by object *pdbs*. With the default argument `rm.gaps=TRUE` unaligned atoms 
#' are omitted from output:

#+ example1g0, cache=TRUE
modes <- nma.pdbs(pdbs, rm.gaps=TRUE, ncore=4)

#'
#' The *modes* object of class *enma* contains aligned normal mode data including fluctuations,
#' RMSIP data, and aligned eigenvectors. A short summary of the *modes* object can be obtain by
#' calling the function **print()**, and the aligned fluctuations can be plotted with function
#' **plot.enma()**. 

#+ example1g1, eval=TRUE,
print(modes)

#+ example1g2, eval=FALSE,
# plot modes fluctuations
plot(modes, pdbs=pdbs, col=grps)

#' **plot.enma()** also facilitates the calculation and visualization of sequence conservation
#' (use argument `conservation=TRUE`), and fluctuation variance (use argument `variance=TRUE`).
#' Note that `conservation` can also be a numeric vector, e.g. manually calculated sequence conservation: 


#+ example1g3, eval=TRUE,
#+ fig1-5, fig.cap="Flexibility profiles and sequence conservation of bacterial DHFR. Upper panel shows the modes fluctuations (calculated with **nma.pdbs()**) of the 81 collected DHFR structures. Line colors depict the different species (E.coli: black; M.profunda: red; B. anthracis green; Y. pestis blue). Lower panel show the sequence conservation calculated based on the results from blasting the swissprot database. Note the low sequence conservation for the two major loops showing large flexibility in the middle of the sequence."

# mode fluctuations and sequence conservation from the swissprot blast
plot(modes, pdbs=pdbs, conservation=cons[gaps.aln$f.inds], col=grps)

#'
#' ### Group by similarity of normal modes
#' The similarity of structural dynamics is calculated by RMSIP based on the 10 lowest
#' frequency normal modes. The *rmsip* values are pre-calculated in the *modes* object
#' and can be accessed through the attribute `modes$rmsip`. The rmsip matrix facilitates
#' clustering of structures with similar flexibility pattern:

hc.nma <- hclust(as.dist(1 - modes$rmsip))
grps.nma <- cutree(hc.nma, k=4)

#+ fig1-6, fig.cap="RMSIP-based clustering of DHFR structures. The heatmap shows the root mean square inner product (RMSIP) for all structure pairs. The structures can be dividen into groups based on their pairwise RMSIP values (row-side colors)."
heatmap((1 - modes$rmsip), labCol = ids, labRow = mynames, 
        symm = TRUE, distfun=as.dist,
        RowSideColors=as.character(grps.nma),
        ColSideColors=as.character(grps),
        col=bwr.colors(10))

#'
#' A function call to **rmsip()** facilitates a more detailed comparison of two modes subspaces as it provides
#' the overlap (squared dot product) of all pairwise mode vectors. 

# Compare the modes of two B.anthracis structures members showing low RMSIP
r1 <- rmsip(modes$U.subspace[,,which(grps.nma==3)[5]],
            modes$U.subspace[,,which(grps.nma==4)[1]])

# Compare the modes of E.coli and B.anthracis structure
r2 <- rmsip(modes$U.subspace[,,which(grps.nma==3)[1]],
            modes$U.subspace[,,which(grps.nma==1)[1]])


#+ fig1-7a, fig.cap="Overlap and RMSIP of (A) two *B.anthracis* structures, and (B) an *E.coli* and a  *B.anthracis* structure.", fig.width=8, fig.height=4.5,
par(mfrow=c(1,2))
plot(r1, xlab=mynames[which(grps.nma==3)[5]], ylab=mynames[which(grps.nma==4)[1]] )
plot(r2, xlab=ids[which(grps.nma==3)[1]], ylab=ids[which(grps.nma==1)[1]] )


#'
#' ### Fluctuation analysis
#' Comparing the mode fluctuations of two groups of structures 
#' can reveal specific regions of distinc flexibility patterns. Below we focus on the
#' *E.coli* and *B.anthracis* structures in our dataset, and use argument `signif=TRUE`
#' to determine significant differences of fluctuations between the two groups:

#+ example1i, cache=TRUE,
# show only grps 1 and 3
cols <- grps
cols[which(grps!=1 & grps!=3)]=NA

#+ fig1-8, fig.cap="Comparison of mode fluctuations between *E.coli* (black), and  *B.anthracis* (green). Residues showing significant differences between the two groups are marked with blue regions.", fig.height=4.5,
ind <- plot(modes, pdbs=pdbs, col=cols, signif=TRUE)

#' Map the differences in fluctuation profiles onto the protein structure:
gaps.res <- gap.inspect(pdbs$ali)
gaps.pos <- gap.inspect(pdbs$xyz)

chain <- rep("A", length(gaps.res$f.inds))
chain[ind$signif]="B"
write.pdb(xyz=pdbs$xyz[1,gaps.pos$f.inds], b=modes$fluctuations[1,],
          chain=chain,
          resno=pdbs$resno[1, gaps.res$f.inds],
          resid=pdbs$resno[1, gaps.res$f.inds])
         

#'
#' ### Visualize modes
#' A function call to **mktrj.enma()** will generate a trajectory pdb file for the visualization of
#' a specific normal mode for one of the structures in the *pdbs* object. This allows for a visual
#' comparison of the calculated normal modes. In this particular example the similarity of the
#' collective modes are striking: 

#+ example1j, cache=TRUE, results='hide', eval=FALSE,
mktrj.enma(modes, pdbs, mode=1, which(grps==1)[1]) ## E coli
mktrj.enma(modes, pdbs, mode=1, which(grps==2)[1]) ## M. profunda
mktrj.enma(modes, pdbs, mode=1, which(grps==3)[1]) ## B. anthracis
mktrj.enma(modes, pdbs, mode=1, which(grps==4)[1]) ## Y. pestis


#' ![Comparison of modes fluctuations profiles of *E.coli* (left panel) and *B.anthracis* (rigth panel).](figure/modes-compare.png)


#'
#' Another useful functionality is the identification of dynamic domains - i.e. groups of residues
#' moving as coherent units. A function call to **geostas()** will attempt to identify such groups
#' based on trajectory data. We therefore make a trajectory consisting of the first 5 modes
#' by interpolating along the mode vectors. The resulting domain assignment are stored in the `gs$grps`
#' attribute, and can be visualized by writing to a PDB file:


#+ example1j2, cache=TRUE, results='hide',
# Find domains for E.coli DHFR
struct <- which(grps==1)[1]
trj <- NULL
for(i in 1:5)
  trj=rbind(trj, mktrj.enma(modes, pdbs, mode=i, struct))

gs <- geostas(trj)
write.pdb(pdb=NULL, xyz=trj, chain=gs$grps, file="R1.pdb")

# Find domains for B.anthracis DHFR
struct <- which(grps==3)[1]
trj <- NULL
for(i in 1:5)
  trj=rbind(trj, mktrj.enma(modes, pdbs, mode=i, struct))

gs <- geostas(trj)
write.pdb(pdb=NULL, xyz=trj, chain=gs$grps, file="R3.pdb")



#' ### RMSD-based clustering
#' As a supplement we show that grouping by pairwise RMSD values provides slightly different results
#' than grouping by RMSIP from NMA. In h

#+ example1g3, eval=TRUE, warning=FALSE,
rd <- rmsd(pdbs)
hc.rd <- hclust(as.dist(rd))
grps.rd <- cutree(hc.rd, k=4)

#+ fig1-7b, fig.cap="RMSD-based clustering reveals slightly different grouping than the RMSIP-based clustering."
heatmap(rd, labCol = ids, labRow = mynames, 
        symm = TRUE, distfun=as.dist,
        RowSideColors=as.character(grps.nma),
        ColSideColors=as.character(grps.rd),
        col=bwr.colors(10))


#'
#' ### Cross-correlation analysis

#+ example1-dccm, cache=TRUE,
##modes <- nma.pdbs(pdbs, ncore=4, rm.gaps=FALSE)
cij <- dccm(modes, ncore=4, na.rm=FALSE)

#+ fig1-9, fig.cap="E.coli DCCM.", fig.height=6, fig.width=6, 
plot.dccm(cij$all.dccm[,,1], contour=F, col.regions=bwr.colors(200), at=seq(-1,1,by=0.01) )

#+ fig1-10, fig.cap="B.anthracis DCCM.", fig.height=6, fig.width=6, 
plot.dccm(cij$all.dccm[,,6], contour=F, col.regions=bwr.colors(200), at=seq(-1,1,by=0.01) )




#'
#' ### Principal Component Analysis 
#' Finally, a principal component analysis (PCA) can be carried out to investigate the major
#' conformational variability of the collected ensemble of structures. User function **pca.xyz()**
#' with `pdbs$xyz` as argument (note that `pdbs$xyz` is already superimposed from the functioncall
#' **pdbaln()** above). 

#+ example1-pca-1, warning=FALSE, cache=TRUE, 
gaps.res <- gap.inspect(pdbs$ali)
gaps.pos <- gap.inspect(pdbs$xyz)

pc.xray <- pca.xyz(pdbs$xyz[, gaps.pos$f.inds])

#+ fig1-11, fig.cap="Principal component analysis of DHFR X-ray structures. The plot shows the relationship between differrent conforms in terms of their major structural displacements (i.e. along the the first and second principal components). Each dot (structure) is colored according to its specie."
plot(pc.xray$z[,1:2], pch=16, xlab="PC 1", ylab="PC 2", col=grps)

# visualize the first principal component
a <- mktrj.pca(pc.xray, pc=1, file="pc1.pdb",
               resno = pdbs$resno[1, gaps.res$f.inds],
               resid = aa123(pdbs$ali[1, gaps.res$f.inds]) )

# compare the PCA and NMA
rmsip(pc.xray$U, modes$U.subspace[,,5])


#'
#' ### Compare with MD simulation
#+ example1-pca-2, warning=FALSE, cache=TRUE, results='hide',
mddir <- "01_md-simulations/01_1rx2"
pdbmd <- paste(mddir, "01_build", "prot.pdb", sep="/")

pdb <- read.pdb(pdbmd)
ca.inds <- atom.select(pdb, 'calpha')
pdb <- trim.pdb(pdb, ca.inds)

trj <- read.ncdf(paste(mddir, "03_prod", "prod_noWAT.nc", sep="/"), at.sel=ca.inds)

md.inds <- add.pdb(pdb, pdbs)
xyz=fit.xyz(fixed=pdbs$xyz[1, gaps.pos$f.inds],
  mobile=trj[, md.inds$xyz])

proj <- pca.project(xyz, pc.xray)
cols <- densCols(proj[,1:2])

#+ fig1-12, fig.cap="Projection of MD conformers onto the X-ray PC space."
plot(proj[,1:2], col=cols, pch=16,
     xlim=range(pc.xray$z[,1]), ylim=range(pc.xray$z[,2]))
points(pc.xray$z[,1:2], col=1, pch=1, cex=1.1)
points(pc.xray$z[,1:2], col=grps, pch=16)

# PCA of the MD trajectory
pc.md <- pca.xyz(xyz)

# compare MD-PCA and NMA
rmsip(pc.md$U, modes$U.subspace[,,5])





#' ## Document Details
#' This document is shipped with the Bio3D package in both R and PDF formats. All code can be extracted and automatically executed to generate Figures and/or the PDF with the following commands:

#+ close, include=TRUE, eval=FALSE
library(knitr)
spin('Bio3D_nma-dhfr-partI.r')
system("pandoc -o Bio3D_nma-dhfr-partI.pdf Bio3D_nma-dhfr-partI.md")

#'
#' ## Information About the Current Bio3D Session
#'
print(sessionInfo(), FALSE)
