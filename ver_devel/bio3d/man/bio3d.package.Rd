\name{bio3d-package}
\alias{bio3d-package}
\alias{bio3d}
\docType{package}
\title{
Biological Structure Analysis
}
\description{
Utilities for the analysis of protein structure and sequence data. 
}
\details{
\tabular{ll}{
Package: \tab bio3d\cr
Type: \tab Package\cr
Version: \tab 2.0\cr
Date: 2013-11-08\cr
License: \tab GPL version 2 or newer\cr
URL: \tab \url{http://thegrantlab.org/bio3d/}\cr
}

Features include the ability to read and write structure
(\code{\link{read.pdb}},  \code{\link{write.pdb}}, 
\code{\link{read.fasta.pdb}}),  sequence (\code{\link{read.fasta}}, 
\code{\link{write.fasta}}) and dynamics trajectory data
(\code{\link{read.dcd}}, \code{\link{read.ncdf}}, \code{\link{write.ncdf}}).

Perform sequence database searches (\code{\link{blast.pdb}}), atom
summaries (\code{\link{summary.pdb}}),  atom selection
(\code{\link{atom.select}}),  re-orientation (\code{\link{orient.pdb}}),
superposition (\code{\link{rot.lsq}},  \code{\link{fit.xyz}}), rigid core
identification (\code{\link{core.find}},  \code{\link{plot.core}}, 
\code{\link{fit.xyz}}), torsion/dihedral analysis
(\code{\link{torsion.pdb}},  \code{\link{torsion.xyz}}), clustering (via
\code{\link{hclust}}), principal component analysis
(\code{\link{pca.xyz}},  \code{\link{pca.tor}},  \code{\link{plot.pca}},
\code{\link{plot.pca.loadings}},  \code{\link{mktrj.pca}}) and dynamical
cross-correlation analysis (\code{\link{dccm}}, \code{\link{plot.dccm}})
of structure data.

Perform conservation analysis of sequence (\code{\link{seqaln}}, \code{\link{conserv}}, 
\code{\link{seqidentity}},  \code{\link{entropy}},  \code{\link{consensus}})
and structural (\code{\link{pdbaln}}, \code{\link{rmsd}},
\code{\link{rmsf}}, \code{\link{core.find}}) data.

Perform normal mode analysis (elastic network model) with
(\code{\link{nma}}) and (\code{\link{build.hessian}}), ensemble normal
mode analysis (\code{\link{nma.pdbs}}), mode comparison
(\code{\link{rmsip}}) and (\code{\link{overlap}}), atomic fluctuation
prediction (\code{\link{fluct.nma}}), cross-correlation analysis
(\code{\link{dccm.nma}}), cross-correlation visualization (\code{\link{view.dccm}}), 
deformation analysis (\code{\link{deformation.nma}}), and mode visualization
(\code{\link{view.modes}}), (\code{\link{mktrj.nma}}). 

In addition, various utility functions are provided to facilitate
manipulation and analysis of biological sequence and structural data
(e.g. \code{\link{get.pdb}}, \code{\link{get.seq}}, \code{\link{aa123}},
\code{\link{aa321}}, \code{\link{pdbseq}}, \code{\link{aln2html}}, \code{\link{atom.select}},
\code{\link{rot.lsq}},  \code{\link{fit.xyz}}, \code{\link{is.gap}}, \code{\link{gap.inspect}},
\code{\link{orient.pdb}} and \code{\link{pairwise}}).
}
\note{
The latest version and further documentation can be obtained from the
bio3d website:\cr
\url{http://thegrantlab.org/bio3d/}.\cr
\url{http://thegrantlab.org/bio3d/html/}.\cr
\url{http://bitbucket.org/Grantlab/bio3d}.
}
\author{
Barry Grant <bjgrant@umich.edu>
Xin-Qiu Yao <xinqyao@umich.edu>
Lars Skjaerven <larsss@gmail.com>
}
\references{
Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\keyword{ documentation }
\examples{
help(package="bio3d")     # list the functions within the package
#lbio3d()                  # list bio3d function names only

## Or visit:
##   http://thegrantlab.org/bio3d/html/

## See the individual functions for further documentation and examples, e.g.
help(read.pdb)

## Or online:
##    http://thegrantlab.org/bio3d/html/read.pdb.html

\dontrun{
##-- Read a PDB file
pdb <- read.pdb("1BG2")

##
##-- Distance matrix
##
k <- dm(pdb, selection="calpha")
plot(k)

## Extract SEQRES PDB sequence
s <- aa321(pdb$seqres)

## Extract ATOM PDB sequence
s2 <- pdbseq(pdb)

## write a FASTA format sequence file
write.fasta(seqs=seqbind(s, s2), id=c("seqres","atom"), file="eg.fa")

##----------------------------------------##

##
##-- Select alpha carbon atom subset
##
ca.inds <- atom.select(pdb, "calpha")

## Plot of B-factor values along with secondary structure
plot.bio3d(pdb$atom[ca.inds$atom, "b"], sse=pdb, ylab="B-factor")

## Secondary structure assignment with DSSP
#sse <- dssp(pdb)

##----------------------------------------##

##
##-- Torsion angle analysis and basic Ramachandran plot
##
tor <- torsion.pdb(pdb)
plot(tor$phi, tor$psi)

##----------------------------------------##

##
##-- Search for related structures in the PDB database
blast <- blast.pdb( pdbseq(pdb) )
hits  <- plot.blast(blast)
head(hits$hits)

## Download these with function "get.pdb()"
#rawpdbs <- get.pdb( hits$pdb.id, "PDB_downloads")

## Split by chain with function "pdbsplit()"
#pdbsplit(rawpdbs, path="PDB_downloads/split_chain")

## and then align with "pdbaln()" and superpose with "pdbfit()"
#hitfiles <- paste0("PDB_downloads/split_chain", hits$pdb.id, ".pdb")
#pdbs <- pdbaln(hitfiles)
#xyz  <- pdbfit(pdbs)

##----------------------------------------##

##
##-- Read an example FASTA sequence alignment from PFAM
##
infile <- "http://pfam.sanger.ac.uk/family/PF00071/alignment/seed/format?format=fasta"

aln <- read.fasta( infile )

## Entropy and similarity scores for alignment positions
h <- entropy(aln)   
s <- conserv(aln)   # see other"conserv()" options
plot(h$H.norm, typ="h", ylab="Normalized entropy score", col="gray")
points( s, typ="h", col="red")

## Alignment consensus sequence
con <- consensus(aln)
con$seq

## add consensus sequence to conservation plot
ind <- which(s > 0.6)
text(ind, s[ind], labels=con$seq[ind])

## Render the alignment as coloured HTML
aln2html(aln, append=FALSE, file="eg.html")

##----------------------------------------##

##
##-- Read an alignment of sequences and their corresponding structures
##
aln  <- read.fasta( system.file("examples/kif1a.fa", package="bio3d") )
pdbs <- read.fasta.pdb( aln )

##-- DDM: Difference Distance Matrix
a <- dm(pdbs$xyz[2,])
b <- dm(pdbs$xyz[3,])
ddm <- a - b
plot(ddm,key=FALSE, grid=FALSE)


##-- Superpose structures on non gap positions
xyz <- pdbfit(pdbs)

##-- RMSD of non gap positions
gaps <- gap.inspect(pdbs$xyz)
rmsd(pdbs$xyz[, gaps$f.inds])
rmsd(xyz[, gaps$f.inds])


##-- Rigid 'core' identification
core <- core.find(pdbs)
#plot(core)

## Core fit the structures (superpose on rigid zones)
xyz2 <- pdbfit(pdbs, inds=core$c0.5A.xyz)

## Note larger overall RMSD but lower core-residue RMSF
rmsd(xyz2[, gaps$f.inds])

plot(rmsf(xyz), typ="l", col="blue", ylab="RMSF")
points(rmsf(xyz2), typ="l", col="red")



##
##-- PCA of experimental structures
##
# Ignore gap containing positions
gaps.res <- gap.inspect(pdbs$ali)
gaps.pos <- gap.inspect(pdbs$xyz)


##-- Do PCA
pc.xray <- pca.xyz(xyz[, gaps.pos$f.inds])

## Plot results
plot(pc.xray)

plot.pca.loadings(pc.xray$au)

## Write a PC trajectory (for viewing as tube in VMD)
rn <- pdbs$resno[1, gaps.res$f.inds]
rd <- aa123(pdbs$ali[1, gaps.res$f.inds])

#a <- mktrj.pca(pc.xray, pc=1, resno =rn, resid = rd, file="pc1.pdb")
#b <- mktrj.pca(pc.xray, pc=2, resno =rn, resid = rd, file="pc2.pdb")
c <- mktrj.pca(pc.xray, pc=3, resno =rn, resid = rd, file="pc3.pdb")


##----------------------------------------##

##
##-- Read a CHARMM/X-PLOR/NAMD trajectory file
##
trtfile <- system.file("examples/hivp.dcd", package="bio3d")
trj <- read.dcd(trtfile)

## Read the starting PDB file to determine atom correspondence
pdbfile <- system.file("examples/hivp.pdb", package="bio3d")
pdb <- read.pdb(pdbfile)

## Fit trj on PDB based on residues 23 to 31 and 84 to 87 in both chains
##inds <- atom.select(pdb,"///23:31,84:87///CA/")
inds <- atom.select(pdb, resno=c(23:31,84:87), elety="CA")
fit.xyz <- fit.xyz(pdb$xyz, trj, fixed.inds=inds$xyz, mobile.inds=inds$xyz)

##-- RMSD of trj frames from PDB
r <- rmsd(a=pdb, b=fit.xyz)

##-- PCA of trj
pc.trj <- pca.xyz(fit.xyz)

## Plot PCA results
plot(pc.trj)

## Examine residue-wise contributions to PCs
plot.pca.loadings(pc.trj$au)

## cluster in PC1 subspace
hc <- hclust(dist(pc.trj$z[,1]))
plot(pc.trj, col=cutree(hc, k=2))


## Write PC trajectory for viewing as tube in VMD
a <- mktrj.pca(pc.trj, pc=1, file="pc1_trj.pdb")

}
## other examples include: normal mode analysis, alignment, clustering etc...

}

