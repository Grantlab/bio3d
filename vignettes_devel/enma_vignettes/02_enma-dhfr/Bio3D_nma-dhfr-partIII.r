#' # Supplementary Data: Enhanced Methods for Normal Mode Analysis with Bio3D
#' ## Lars Skjaerven, Xin-Qiu Yao & Barry J. Grant

#+ setup, include=FALSE
opts_chunk$set(dev='pdf')

#+ preamble, include=FALSE, eval=FALSE
library(knitr)
spin('Bio3D_nma-dhfr-partIII.r')
system("pandoc -o Bio3D_nma-dhfr-partIII.pdf Bio3D_nma-dhfr-partIII.md")

#' ## Background:
#' Bio3D is an R package that provides interactive tools for structural bioinformatics. 
#' The primary focus of Bio3D is the analysis of bimolecular structure, sequence and simulation data.

#'
#' ## Part III:  Ensemble NMA across multiple species of DHFR
#' In this vignette we extend the analysis from Part I by including 
#' a more extensive search of distant homologues within the 
#' DHFR family. Based on a HMMER search we identify and collect 
#' protein species down to a pairwise sequence identity of 19\%. 
#' Normal modes analysis (NMA) across these species reveals  
#' a remarkable similarity of the fluctuation profiles, but also 
#' features which are characteristic to specific species.



#+ example3, results="hide", cache=TRUE,
library(devtools)
load_all("~/workspace/bio3d/ver_devel/bio3d")


#' ### HMMER search
#' Below we use the sequence of *E.coli* DHFR to perform an initial search against
#' the Pfam HMM database with function **hmmer()**. The arguments `type` and `db`
#' specifies the type of hmmer search and the database to search, respectively. In this particular
#' example, our query sequence is searched against the Pfam profile HMM library
#' (arguments `type=hmmscan` and `db=pfam`) to identify its respecitve protein family. 
#' The to **hmmer()** will return a data frame object containing the Pfam accession ID (`$acc`), description of the 
#' identified family (`$desc`), family name (`$name`), etc. 

#+ example3a, results="hide", cache=TRUE, warning=FALSE, 
# get sequence of Ecoli DHFR
seq <- get.seq("1rx2_A")

#+ example3a2, cache=TRUE,
# search the Pfam database
pfam <- hmmer(seq, type="hmmscan", db="pfam")
pfam$name
pfam$desc

#' 
#' **Sidenote:** The **hmmer()** function facilitates four different types of searches at a multitude of databases. 
#' Use function `help(hmmer)` for a complete overview of the different options. 

#' 
#' Having identified the Pfam entry of our query protein we can use function **pfam()** to fetch the curated sequence alignment of
#' the DHFR family. Use function **print.fasta()** to print a short summary of the downloaded sequence alignment to the screen.
#' Note that if argument `alignment=TRUE` the sequence alignmenment itself will be written to screen. 

#+ example3a3, cache=TRUE,
# download pfam alignment for family
pfam.aln <- pfam(pfam$acc[1])
print(pfam.aln, alignment=FALSE)

#'
#' The next hmmer search builds a profile HMM from the Pfam multiple sequence alignment
#' and uses this HMM to search against a target sequence database (use `type=hmmsearch`).
#' In this case our target sequence database is the PDB (`db=pdb`), but there are also other options 
#' such as *Swissprot* and *UniProt*. 


#+ example3a4, results="hide", cache=TRUE, warning=FALSE, 
# use Pfam alignment in search 
hmm <- hmmer(pfam.aln, type="hmmsearch", db="pdb")

#'
#' Function **plot.hmmer()** (the equivalent to **plot.blast()**) provides a quick overview of the
#' search results, and can aid in the identification of a sensible hit similarity threshold.
#' The normalized scores (-log(E-Value)) are shown in the upper panel, and the lower panel provides an overview of
#' the kingdom and specie each hit are associated with. Here we specify a cutoff of 55 yielding 104 hits:


#+ fig3-2, fig.cap="Overview  of hits obtained from the HMMER search. Upper panel shows the normalized scores. Lower panel the scores and hits are colored according to their respective kingdom (background colors) and specie (foreground barplot).", fig.height=5,
hits <- plot.hmmer(hmm, cutoff=90) ##55

ids <- hits$acc
species <- hmm$species[hits$inds]

#+ example3a5, cache=TRUE, warning=FALSE, 
# check out what species we got
print(unique(species))


#'
#' Having identified relevant PDB structures through the hmmer search
#' we proceed by fetching and pre-processing
#' the PDB files with functions **get.pdb()** and **pdbsplit()**.

#' 
#' As in the previous vignette,
#' we are interested in protein structures without missing in-structure residues,
#' and we also want to limit the number of identifical conformers:

#+ example3b, results="hide", cache=TRUE, warning=FALSE,
# fetch and split PDBs
raw.files <- get.pdb(ids, path = "raw_pdbs")
files <- pdbsplit(raw.files, ids = ids, path = "raw_pdbs/split_chain",
                  het2atom=FALSE, ncore=4)
pdbs <- pdbaln(files)

# exclude a DHFR-TS fusion protein
inds <- grep("1qzf", pdbs$id, invert=TRUE)
pdbs <- pdbs.filter(pdbs, row.inds=inds)

# exclude structures with missing residues
conn <- inspect.connectivity(pdbs, cut=4.05)
pdbs <- pdbs.filter(pdbs, row.inds=which(conn))
          
# exclude conformational redundant structures
rd <- rmsd.filter(pdbs$xyz, cutoff=0.1, fit=TRUE)
pdbs <- pdbs.filter(pdbs, row.inds=rd$ind)

#' 
#' In this particular case
#' a standard sequence alignment (e.g. through function **pdbaln()** or **seqaln()**) is not sufficient for a proper alignment.
#' We will therefore make use of the Pfam profile alignment, and align our selected PDBs to this
#' using argument `profile` to function **seqaln()**. Subsequently, we re-read the fasta file, and use function
#' **read.fasta.pdb()** to obtain aligned C-alpha atom data (including coordinates etc.) for the PDB ensemble:

#+ example3c, results="hide", cache=TRUE, warning=FALSE,
# align pdbs to Pfam-profile alignment
aln <- seqaln(pdbs, profile=pfam.aln, exefile="clustalo", extra.args="--dealign")

# store only PDBs in alignment
aln$ali=aln$ali[1:length(pdbs$id),]
aln$id=aln$id[1:length(pdbs$id)]

# re-read PDBs to match the new alignment
pdbs <- read.fasta.pdb(aln)

# exclude gap-only columns
pdbs=pdbs.filter(pdbs)

# refit coordinates
pdbs$xyz <- pdbfit(pdbs, outpath="raw_pdbs/flsq/")

# fetch IDs again
ids = unlist(strsplit(basename(pdbs$id), split=".pdb"))
species = hmm$species[hmm$acc %in% ids]

#+ example3d, results="hide", cache=TRUE, warning=FALSE,
# labels for annotating plots
mynames <- paste(substr(species, 1,1), ". ", 
                 lapply(strsplit(species, " "), function(x) x[2]), sep="")
mynames

#'
#' The *pdbs* object now contains *aligned* C-alpha atom data, including Cartesian coordinates,
#' residue numbers, residue types, and B-factors. The sequence alignment is also stored by default
#' to the FASTA format file 'aln.fa' (to view this you can use an alignment viewer such as SEAVIEW,
#' see _Requirements_ section above). 


#'
#' ### Sequence conservation analysis
#' Function **seqidentity()** can be used to calculate the sequence identity for the PDBs ensemble. 
#' Below we also print a summary of the calculated sequence identities, and perform a clustering
#' of the structures based on  sequence identity:

#+ example3e, warning=FALSE, cache=TRUE, 
seqide <- seqidentity(pdbs)
summary(c(seqide))
hc <- hclust(as.dist(1-seqide))
grps.seq <- cutree(hc, k=5)

#+ fig3-3, fig.cap="Sequence-based clustering of the collected structures.", fig.height=5,
plot(hc, hang=-1, labels=species, cex=0.5)


#'
#' ### Normal modes analysis
#' Function **nma.pdbs()** will calculate the normal modes of each protein structures stored
#' in the *pdbs* object. The normal modes are calculated on the full structures as provided
#' by object *pdbs*. Use argument `rm.gaps=FALSE` to visualize fluctuations also of un-aligned
#' residues:

#+ example3f, cache=TRUE
modes <- nma.pdbs(pdbs, rm.gaps=FALSE, ncore=4)

#'
#' The *modes* object of class *enma* contains aligned normal mode data including fluctuations,
#' RMSIP data (only when`rm.gaps=FALSE`), and aligned eigenvectors. A short summary of the *modes* object can be obtain by
#' calling the function **print()**, and the aligned fluctuations can be plotted with function
#' **plot.enma()**:
print(modes)

#+ fig3-4, fig.cap="Flexibility profiles and sequence conservation. The figure shows the modes fluctuations colored according their sequence identity. The lower panel shows the sequence conservation for the PDBs."
plot(modes, pdbs=pdbs, ylim=c(0,2), col=grps.seq, conservation=TRUE)


#'
#' ### Visualize modes
#' A function call to **mktrj.enma()** will generate a trajectory PDB file for the visualization of
#' a specific normal mode for one of the structures in the *pdbs* object. This allows for a visual
#' comparison of the calculated normal modes. Below we make a PDB trajectory of the first mode (argument `mode=1`)
#' of 4 relevant species (e.g. argument `ind=1`). Note that we use **grep()** to fetch the
#' indices (in the *modes* and *pdbs* objects) of the relevant species:

#+ example3h, cache=TRUE, results='hide', eval=FALSE,
# E.coli
mktrj.enma(modes, pdbs, mode=1, ind=grep("coli", species)[1])

# T. mycobacterium
mktrj.enma(modes, pdbs, mode=1, ind=grep("tuberculosis", species)[1])

# C. albicans
mktrj.enma(modes, pdbs, mode=1, ind=grep("albicans", species)[1])

# H. sapiens
mktrj.enma(modes, pdbs, mode=1, ind=grep("sapiens", species)[1]) 


#' ![Mode comparison of (A) *E.coli*, (B) *Mycobacterium tuberculosis*, (C) *C.albicans*, and (D) *H.sapiens*. The trajectories are made with function **mktrj.enma** and visualized in PyMol.](visualize-4hfrs.png)


#'
#' ### Group by similarity of normal modes
#' The similarity of structural dynamics is calculated by *RMSIP* based on the 10 lowest
#' frequency normal modes. The RMSIP values are pre-calculated in the *modes* object (when `rm.gaps=TRUE`)
#' and can be accessed through the attribute `modes$rmsip`. The matrix of pairwise RMSIP values facilitates
#' clustering of structures based on their flexibility pattern:

#+ example3i, cache=TRUE,
modes <- nma.pdbs(pdbs, ncore=4, rm.gaps=TRUE)

# RMSIP-based clustering
hc.nma <- hclust(as.dist(1 - modes$rmsip))
grps.nma <- cutree(hc.nma, k=5)

#+ fig3-5, fig.cap="RMSIP-based clustering of DHFR structures. The heatmap shows the root mean square inner product (RMSIP) for all structure pairs. The structures can be dividen into groups based on their pairwise RMSIP values (row-side colors), and pairwise sequence identity (column-side colors). "
# RMSIP heatmap
heatmap((1 - modes$rmsip), labCol = ids, labRow = mynames, 
        symm = TRUE, distfun=as.dist,
        RowSideColors=as.character(grps.nma),
        ColSideColors=as.character(grps.seq),
        col=bwr.colors(10))

#'
#' ### Fluctuation analysis
#' **Sidenote:** Note that the gap-regions are excluded from `modes$fluctuations` when `rm.gaps=TRUE`
#' in the modes calculation:

#+ fig3-7, fig.cap="Fluctuation profiles. Note that the gap-regions are excluded from `modes$fluctuations` when `rm.gaps=TRUE` in the modes calculation. ", fig.height=4,
plot(modes, pdbs=pdbs, col=grps.seq)

#'
#' ### PCA
gaps.res <- gap.inspect(pdbs$ali)
gaps.pos <- gap.inspect(pdbs$xyz)

pc.xray <- pca.xyz(pdbs$xyz[, gaps.pos$f.inds])

#+ fig3-9, fig.cap="Principal component analysis of DHFR X-ray structures. The plot shows the relationship between differrent conforms in terms of their major structural displacements (i.e. along the the first and second principal components). Each dot (structure) is colored according to its specie."
plot(pc.xray$z[,1:2], pch=16, xlab="PC 1", ylab="PC 2", col=grps.nma)

# visualize the first principal component
a <- mktrj.pca(pc.xray, pc=1, file="pc1.pdb",
               resno = pdbs$resno[1, gaps.res$f.inds],
               resid = aa123(pdbs$ali[1, gaps.res$f.inds]) )



#' ## Document Details
#' This document is shipped with the Bio3D package in both R and PDF formats. All code can be extracted and automatically executed to generate Figures and/or the PDF with the following commands:

#+ close, include=TRUE, eval=FALSE
library(knitr)
spin('Bio3D_nma-dhfr-partIII.r')
system("pandoc -o Bio3D_nma-dhfr-partIII.pdf Bio3D_nma-dhfr-partIII.md")

#'
#' ## Information About the Current Bio3D Session
#'
print(sessionInfo(), FALSE)



