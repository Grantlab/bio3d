#' # Supporting Material S4
#' # Integrated structural and evolutionary ensemble analysis with Bio3D
#' **Lars Skj\ae rven, Xin-Qiu Yao & Barry J. Grant**

#+ setup, include=FALSE
opts_chunk$set(dev='pdf')

#+ preamble, include=FALSE, eval=FALSE
library(knitr)
spin('Bio3D_nma-dhfr-partII.r')
system("pandoc -o Bio3D_nma-dhfr-partII.pdf Bio3D_nma-dhfr-partII.md")

#' ## Background:
#' Bio3D[^1] is an R package that provides interactive tools for structural bioinformatics. The primary focus of Bio3D is the analysis of bimolecular structure, sequence and simulation data.

#'
#' Normal mode analysis (NMA) is one of the major simulation techniques used to probe large-scale motions in biomolecules. Typical application is for the prediction of functional motions in proteins. Version 2.0 of the Bio3D package now includes extensive NMA facilities (see also the [**NMA Vignette**](http://thegrantlab.org/bio3d/tutorials)). These include a unique collection of multiple elastic network model force-fields, automated ensemble analysis methods, and variance weighted NMA. Here we provide an in-depth demonstration of ensemble NMA with working code that comprise complete executable examples[^2].

#'
#' [^1]: The latest version of the package, full documentation and further vignettes (including detailed installation instructions) can be obtained from the main Bio3D website: [http://thegrantlab.org/bio3d/](http://thegrantlab.org/bio3d/)
#'
#' [^2]: This document contains executable code that generates all figures contained within this document. See help(vignette) within R for full details. 


#'
#' #### Requirements:
#' Detailed instructions for obtaining and installing the Bio3D package on various platforms can be found in the [**Installing Bio3D Vignette**](http://thegrantlab.org/bio3d/download/download.html) available both on-line and from within the Bio3D package. In addition to Bio3D the _MUSCLE_ and _CLUSTALO_ multiple sequence alignment programs (available from the [muscle home page](http://www.drive5.com/muscle/) and [clustalo home page](http://www.clustal.org/omega)) must be installed on your system and in the search path for executables. Please see the installation vignette for further details.


#'
#' ## Part II:  Ensemble NMA across multiple species of DHFR
#' In this vignette we extend the analysis from Part I by including 
#' a more extensive search of distant homologues within the
#' DHFR family. Based on a HMMER search we identify and collect
#' protein species down to a pairwise sequence identity of 21\%.
#' Normal modes analysis (NMA) across these species reveals
#' a remarkable similarity of the fluctuation profiles, but also
#' features which are characteristic to specific species.


#'
#' ### HMMER search distantly related DHFR species
#' Below we use the sequence of *E.coli* DHFR to perform an initial search against
#' the Pfam HMM database with function **hmmer()**. The arguments `type` and `db`
#' specifies the type of hmmer search and the database to search, respectively. In this particular
#' example, our query sequence is searched against the Pfam profile HMM library
#' (arguments `type=hmmscan` and `db=pfam`) to identify its respecitve protein family. 
#' The to **hmmer()** will return a data frame object containing the Pfam accession ID (`$acc`),
#' description of the identified family (`$desc`), family name (`$name`), etc. 

#+ getseq, results="hide", cache=TRUE, warning=FALSE,
# load the bio3d package
library(bio3d)

# get sequence of Ecoli DHFR
seq <- get.seq("1rx2_A")

#+ pfam, cache=TRUE, message=FALSE,
# scan the Pfam database for our sequence
pfam <- hmmer(seq, type="hmmscan", db="pfam")
pfam$name
pfam$desc

#' 
#' **Sidenote:** The **hmmer()** function facilitates four different types of searches at a
#' multitude of databases. Use function `help(hmmer)` for a complete overview of the different options. 

#' 
#' Having identified the Pfam entry of our query protein we can use function **pfam()** to fetch the
#' curated sequence alignment of the DHFR family. Use function **print.fasta()** to print a short summary
#' of the downloaded sequence alignment to the screen. Note that if argument `alignment=TRUE`
#' the sequence alignmenment itself will be written to screen. 

#+ pfam2, cache=TRUE,
# download pfam alignment for the DHFR family
pfam.aln <- pfam(pfam$acc[1])
print(pfam.aln, alignment=FALSE)

#'
#' The next hmmer search builds a profile HMM from the Pfam multiple sequence alignment
#' and uses this HMM to search against a target sequence database (use `type=hmmsearch`).
#' In this case our target sequence database is the PDB (`db=pdb`), but there are also other options 
#' such as *Swissprot* and *UniProt*. 


#+ hmmer, results="hide", cache=TRUE, warning=FALSE, 
# use Pfam alignment in search 
hmm <- hmmer(pfam.aln, type="hmmsearch", db="pdb")

#'
#' Function **plot.hmmer()** (the equivalent to **plot.blast()**) provides a quick overview of the
#' search results, and can aid in the identification of a sensible hit similarity threshold.
#' The normalized scores (-log(E-Value)) are shown in the upper panel, and the lower panel provides an
#' overview of the kingdom and specie each hit are associated with. Here we specify a cutoff
#' of 90 yielding 517 hits:


#+ fig3-2, fig.cap="Overview  of hits obtained from the HMMER search. Upper panel shows the normalized scores. Lower panel the scores and hits are colored according to their respective kingdom (background colors) and specie (foreground barplot).", fig.height=5,
hits <- plot.hmmer(hmm, cutoff=90)

ids <- hits$acc
species <- hmm$species[hits$inds]

#+ species, cache=TRUE, warning=FALSE, 
# collected species
print(unique(species))


#'
#' ### Retrieve and process structures from the PDB
#' Having identified relevant PDB structures through the hmmer search
#' we proceed by fetching and pre-processing
#' the PDB files with functions **get.pdb()** and **pdbsplit()**.

#' 
#' As in the previous vignette,
#' we are interested in protein structures without missing in-structure residues,
#' and we also want to limit the number of identifical conformers:

#+ pdbs, results="hide", cache=TRUE, warning=FALSE, message=FALSE,
# fetch and split PDBs
raw.files <- get.pdb(ids, path = "raw_pdbs", gzip=TRUE)
files <- pdbsplit(raw.files, ids = ids,
                  path = "raw_pdbs/split_chain", ncore=4)
pdbs.all <- pdbaln(files)

# exclude DHFR-TS fusion protein
excl.inds <- unlist(lapply(c("1qzf", "1sej"), grep, pdbs.all$id))
pdbs <- pdbs.filter(pdbs.all, row.inds=-excl.inds)

# exclude structures with missing residues
conn <- inspect.connectivity(pdbs, cut=4.05)
pdbs <- pdbs.filter(pdbs, row.inds=which(conn))
          
# exclude conformational redundant structures
rd <- rmsd.filter(pdbs$xyz, cutoff=0.25, fit=TRUE, ncore=4)
pdbs <- pdbs.filter(pdbs, row.inds=rd$ind)

#' 
#' In this particular case a standard sequence alignment (e.g. through function **pdbaln()**
#' or **seqaln()**) is not sufficient for a correct alignment.
#' We will therefore make use of the Pfam profile alignment, and align our selected PDBs to this
#' using argument `profile` to function **seqaln()**. Subsequently, we re-read the fasta file, and use function
#' **read.fasta.pdb()** to obtain aligned C-alpha atom data (including coordinates etc.) for the PDB ensemble:

#+ realign, results="hide", cache=TRUE, warning=FALSE,
# align pdbs to Pfam-profile alignment
aln <- seqaln(pdbs, profile=pfam.aln, exefile="clustalo", extra.args="--dealign")

# store only PDBs in alignment
aln$ali=aln$ali[1:length(pdbs$id),]
aln$id=aln$id[1:length(pdbs$id)]

# re-read PDBs to match the new alignment
pdbs <- read.fasta.pdb(aln)

# exclude gap-only columns
pdbs = pdbs.filter(pdbs)

# refit coordinates
pdbs$xyz = pdbfit(pdbs)

#+ pdbfit2, eval=FALSE
# refit coordinates, and write PDBs to disk
pdbs$xyz <- pdbfit(pdbs, outpath="flsq/")

#+ ids, cache=TRUE, warning=FALSE,
# fetch IDs again
ids = unlist(strsplit(basename(pdbs$id), split="\\.pdb"))
species = hmm$species[hmm$acc %in% ids]


#+ labels, cache=TRUE, warning=FALSE,
# labels for annotating plots
labs <- paste(substr(species, 1,1), ". ", 
                 lapply(strsplit(species, " "), function(x) x[2]), sep="")
print(unique(labs))

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

#+ seqide, warning=FALSE, cache=TRUE, 
seqide <- seqidentity(pdbs)
summary(c(seqide))
hc <- hclust(as.dist(1-seqide))
grps.seq <- cutree(hc, h=0.6)

#+ fig3-3, fig.cap="Clustering of collected structures based on sequence identity.", fig.height=5,
plot(hc, hang=-1, labels=labs, cex=0.25)


#'
#' ### Normal modes analysis
#' Function **nma.pdbs()** will calculate the normal modes of each protein structures stored
#' in the *pdbs* object. The normal modes are calculated on the full structures as provided
#' by object *pdbs*. Use argument `rm.gaps=FALSE` to visualize fluctuations also of un-aligned
#' residues:

#+ modes1, cache=TRUE, message=FALSE
modes <- nma.pdbs(pdbs, rm.gaps=FALSE, ncore=4)

#'
#' The *modes* object of class *enma* contains aligned normal mode data including fluctuations,
#' RMSIP data (only when`rm.gaps=FALSE`), and aligned eigenvectors. A short summary of the *modes* object can be obtain by
#' calling the function **print()**, and the aligned fluctuations can be plotted with function
#' **plot.enma()**:
print(modes)

#+ fig3-4, fig.cap="Flexibility profiles and sequence conservation. The figure shows the modes fluctuations colored according their sequence identity. The lower panel shows the sequence conservation for the PDBs. The plot is generated with function **plot.enma()** with argument `conservation=TRUE`."
plot(modes, pdbs=pdbs, ylim=c(0,2), col=grps.seq, conservation=TRUE)

#'
#' In some cases it can be difficult to interpret the fluctuation plot when all lines are plotted on top of each other. Argument `spread=TRUE` adds a small gap between grouped fluctuation profiles. Use this argument in combination with a new groups (`grps`) variable to function **plot.enma()**:

#+ fig3-5, fig.cap="Flexibility profiles for three selected species (*E.coli* (black), *H.sapiens* (red), and *C.albicans* (green)). The plot is generated with function **plot.enma()** with argument `spread=TRUE`.", fig.height=4.5,

grps <- rep(NA, length(grps.seq))
grps[grep("coli", labs)]=1
grps[grep("sapiens", labs)]=2
grps[grep("albicans", labs)]=3

plot(modes, pdbs=pdbs, col=grps, spread=TRUE)


#'
#' ### Visualize modes
#' A function call to **mktrj.enma()** will generate a trajectory PDB file for the visualization of
#' a specific normal mode for one of the structures in the *pdbs* object. This allows for a visual
#' comparison of the calculated normal modes. Below we make a PDB trajectory of the first mode (argument `mode=1`)
#' of 4 relevant species (e.g. argument `ind=1`). Note that we use **grep()** to fetch the
#' indices (in the *modes* and *pdbs* objects) of the relevant species:

#+ example2h, cache=TRUE, results='hide', eval=FALSE,

inds <- c(grep("coli", species)[1], 
          grep("sapiens", species)[1],
          grep("albicans", species)[1])

# E.coli
mktrj.enma(modes, pdbs, mode=1, ind=inds[1], file="ecoli-mode1.pdb")

# H. sapiens
mktrj.enma(modes, pdbs, mode=1, ind=inds[2], file="hsapiens-mode1.pdb")

# C. albicans
mktrj.enma(modes, pdbs, mode=1, ind=inds[3], file="calbicans-mode1.pdb")


#' ![Mode comparison of *E.coli*, *H.sapiens*, and *C.albicans*. The trajectories are made with function **mktrj.enma** and visualized in PyMol.](figure/visualize-4hfrs.png)


#'
#' ## Document Details
#' This document is shipped with the Bio3D package in both R and PDF formats. All code can be extracted and automatically executed to generate Figures and/or the PDF with the following commands:

#+ close, include=TRUE, eval=FALSE
library(knitr)
spin('Bio3D_nma-dhfr-partII.r')
system("pandoc -o Bio3D_nma-dhfr-partII.pdf Bio3D_nma-dhfr-partII.md")

#'
#' ## Information About the Current Bio3D Session
#'
print(sessionInfo(), FALSE)
