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
Version: \tab 1.1-6\cr
Date: \tab DEVEL VERSION: 2012-11-13\cr
License: \tab GPL version 2 or newer\cr
URL: \tab \url{http://thegrantlab.org/bio3d/}\cr
}

Features include the ability to read and write structure
(\code{\link{read.pdb}},  \code{\link{write.pdb}}, 
\code{\link{read.fasta.pdb}}),  sequence (\code{\link{read.fasta}}, 
\code{\link{write.fasta}}) and dynamics trajectory data
(\code{\link{read.dcd}}, \code{\link{read.ncdf}}, \code{\link{write.ncdf}}).

Perform sequence database searches (\code{\link{blast.pdb}}), atom
summaries (\code{\link{pdb.summary}}),  atom selection
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
\code{\link{identity}},  \code{\link{entropy}},  \code{\link{consensus}})
and structural (\code{\link{pdbaln}}, \code{\link{rmsd}},
\code{\link{rmsf}}, \code{\link{core.find}}) data.

In addition, various utility functions are provided to facilitate
manipulation and analysis of biological sequence and structural data
(e.g. \code{\link{get.pdb}}, \code{\link{get.seq}}, \code{\link{aa123}},
\code{\link{aa321}}, \code{\link{seq.pdb}}, \code{\link{aln2html}}, \code{\link{atom.select}},
\code{\link{rot.lsq}},  \code{\link{fit.xyz}}, \code{\link{is.gap}}, \code{\link{gap.inspect}},
\code{\link{orient.pdb}} and \code{\link{pairwise}}).
}
\note{
The latest version and further documentation can be obtained from the
bio3d website:\cr
\url{http://thegrantlab.org/bio3d/}.\cr
Or its alternate site:\cr
\url{http://mccammon.ucsd.edu/~bgrant/bio3d/}.
}
\author{
Barry Grant <bjgrant@umich.edu>
}
\references{
Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\keyword{ documentation }
\examples{
\dontrun{
help(package="bio3d")  # list the functions within the package
lbio3d()  # list bio3d function names only

## See the individual functions for further documentation and examples.
#help(read.pdb)

##-- Read a PDB file
pdb <- read.pdb("1BG2")

##-- Distance matrix
k <- dm(pdb, selection="calpha")
plot(k)

## Extract SEQRES PDB sequence
s <- aa321(pdb$seqres)

## Extract ATOM PDB sequence
s2 <- seq.pdb(pdb)

## write a FASTA format sequence file
write.fasta(seqs=seqbind(s, s2), id=c("seqres","atom"), file="eg.fa")



## Reorient and write a new PDB file
xyz <- orient.pdb(pdb)
write.pdb(pdb, xyz = xyz, file = "mov1.pdb")


##-- Torsion angle analysis and basic Ramachandran plot
tor <- torsion.pdb(pdb)
plot(tor$phi, tor$psi)


##-- Secondary structure assignment with DSSP
sse <- dssp(pdb, resno=FALSE)

## Plot of B-factor values along with secondary structure
plot.bio3d(pdb$atom[pdb$calpha,"b"], sse=sse, ylab="B-factor")



##-- Search for related structures in the PDB database
blast <- blast.pdb( seq.pdb(pdb) )
hits  <- plot.blast(blast)
head(hits$hits)
## Download these with function "get.pdb()" 
## Split by chain with function "split.pdb()"
## and then align with function "pdbaln()"
rawpdbs <- get.pdb(hits$pdb.id, "PDB_downloads")
split.pdb(rawpdbs)
hitfiles <- paste("split_chain/", hits$pdb.id, ".pdb", sep="")
hitfiles <- hitfiles[-grep("KIN",hitfiles)]
pdbs <- pdbaln(hitfiles)
xyz  <- fit.xyz( pdbs$xyz[1,], pdbs, full.pdbs=TRUE, het=TRUE)

##-- Read an example FASTA alignment file
aln <- read.fasta( system.file("examples/kinesin_xray.fa",package="bio3d") )

## Entropy score for alignment positions
h <- entropy(aln)    # see also the "conserv()" function

# Alignment consensus
con <- consensus(aln)

## Plot of residue conservation (similarity) with secondary structure
plot.bio3d( conserv(aln)[!is.gap(aln$ali[1,])],  sse=sse )

## Render the alignment as coloured HTML
aln2html(aln, append=FALSE, file="eg.html")



##-- Read an alignment of sequences and their corresponding structures
aln  <- read.fasta( system.file("examples/kif1a.fa", package="bio3d") )
pdbs <- read.fasta.pdb( aln, pdbext = ".ent",
               pdb.path = system.file("examples",package="bio3d") )


##-- DDM: Difference Distance Matrix
a <- dm(pdbs$xyz[2,])
b <- dm(pdbs$xyz[3,])
ddm <- a - b
plot(ddm,key=FALSE, grid=FALSE)


##-- Superpose structures on non gap positions
gaps <- gap.inspect(pdbs$xyz)

xyz  <- fit.xyz( fixed  = pdbs$xyz[1,],
                mobile = pdbs$xyz,
                fixed.inds  = gaps$f.inds,
                mobile.inds = gaps$f.inds )


##-- RMSD
rmsd(pdbs$xyz[1,gaps$f.inds], pdbs$xyz[,gaps$f.inds])
rmsd(pdbs$xyz[1,gaps$f.inds], pdbs$xyz[,gaps$f.inds], fit=TRUE)



##-- Rigid 'core' identification
aln <- read.fasta(system.file("examples/kinesin_xray.fa",package="bio3d"))
pdb.path = system.file("examples/",package="bio3d")
pdbs <- read.fasta.pdb(aln, pdb.path = pdb.path, pdbext = ".ent")

core <- core.find(pdbs,
                  #write.pdbs=TRUE,
                  rm.island=TRUE,
                  verbose=TRUE)

## Plot core volume vs length
plot(core)

## Core fit the structures (superpose on rigid zones)
xyz <- fit.xyz( fixed = pdbs$xyz[1,],
                mobile = pdbs,
              #  full.pdbs = TRUE,
              #  pdb.path = pdb.path,
              #  pdbext = ".ent",
              #  outpath = "fitlsq/",
                fixed.inds  = core$c0.5A.xyz,
                mobile.inds = core$c0.5A.xyz)


##-- PCA of structures
cut.seqs <- which(pdbs$id \%in\% c("d1n6mb_","d1ry6a_"))
# Ignore gap containing positions
gaps.res <- gap.inspect(pdbs$ali[-cut.seqs,])
gaps.pos <- gap.inspect(pdbs$xyz[-cut.seqs,])


##-- Do PCA
pc.xray <- pca.xyz(xyz[-cut.seqs, gaps.pos$f.inds])

## Plot results
plot(pc.xray)

#x11()
plot.pca.loadings(pc.xray$U)

## Write PC trajectory
a <- mktrj.pca(pc.xray, pc=1, file="pc1.pdb",
               resno = pdbs$resno[1, gaps.res$f.inds],
               resid = aa123(pdbs$ali[1, gaps.res$f.inds]) )





##-- Read a CHARMM/X-PLOR/NAMD trajectory file
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

plot.pca.loadings(pc.trj)

## Write PC trajectory for viewing as tube in VMD
a <- mktrj.pca(pc.trj, pc=1, file="pc1_trj.pdb")

}
## other examples include: sdm, alignment, clustering etc...

}

