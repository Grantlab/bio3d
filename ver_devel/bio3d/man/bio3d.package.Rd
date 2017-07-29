\name{bio3d-package}
\alias{bio3d-package}
\alias{bio3d}
\docType{package}
\title{ Biological Structure Analysis }
\description{
Utilities for the analysis of protein structure and sequence data.
}
\details{
\tabular{ll}{
Package: \tab bio3d\cr
Type: \tab Package\cr
Version: \tab 2.3-3\cr
Date: \tab 2017-07-29\cr
License: \tab GPL version 2 or newer\cr
URL: \tab \url{http://thegrantlab.org/bio3d/}\cr
}

Features include the ability to read and write structure
(\code{\link{read.pdb}},  \code{\link{write.pdb}},
\code{\link{read.fasta.pdb}}),  sequence (\code{\link{read.fasta}},
\code{\link{write.fasta}}) and dynamics trajectory data
(\code{\link{read.dcd}}, \code{\link{read.ncdf}}, \code{\link{write.ncdf}}).

Perform sequence and structure database searches (\code{\link{blast.pdb}},
\code{\link{hmmer}}), atom summaries (\code{\link{summary.pdb}}),  atom selection
(\code{\link{atom.select}}), alignment (\code{\link{pdbaln}},  \code{\link{seqaln}},
\code{\link{mustang}}) superposition (\code{\link{rot.lsq}},  \code{\link{fit.xyz}}),
\code{\link{pdbfit}}), rigid core identification (\code{\link{core.find}},  \code{\link{plot.core}},
\code{\link{fit.xyz}}), dynamic domain analysis (\code{\link{geostas}}), torsion/dihedral analysis
(\code{\link{torsion.pdb}},  \code{\link{torsion.xyz}}), clustering (via
\code{\link{hclust}}), principal component analysis
(\code{\link{pca.xyz}},  \code{\link{pca.pdbs}}, \code{\link{pca.tor}},  \code{\link{plot.pca}},
\code{\link{plot.pca.loadings}},  \code{\link{mktrj.pca}}), dynamical
cross-correlation analysis (\code{\link{dccm}}, \code{\link{lmi}}, \code{\link{plot.dccm}}) and correlation network analysis (\code{\link{cna}},  \code{\link{plot.cna}}, \code{\link{cnapath}}) of structure data.

Perform conservation analysis of sequence (\code{\link{seqaln}}, \code{\link{conserv}},
\code{\link{seqidentity}},  \code{\link{entropy}},  \code{\link{consensus}})
and structural (\code{\link{pdbaln}}, \code{\link{rmsd}},
\code{\link{rmsf}}, \code{\link{core.find}}) data.

Perform normal mode analysis (\code{\link{nma}}, \code{\link{build.hessian}}), ensemble normal
mode analysis (\code{\link{nma.pdbs}}), mode comparison
(\code{\link{rmsip}}) and (\code{\link{overlap}}), atomic fluctuation
prediction (\code{\link{fluct.nma}}), cross-correlation analysis
(\code{\link{dccm.nma}}), cross-correlation visualization (\code{\link{pymol.dccm}}),
deformation analysis (\code{\link{deformation.nma}}), and mode visualization
(\code{\link{pymol.modes}}, \code{\link{mktrj.nma}}).

In addition, various utility functions are provided to facilitate
manipulation and analysis of biological sequence and structural data
(e.g. \code{\link{get.pdb}}, \code{\link{get.seq}}, \code{\link{aa123}},
\code{\link{aa321}}, \code{\link{pdbseq}}, \code{\link{aln2html}}, \code{\link{atom.select}},
\code{\link{rot.lsq}},  \code{\link{fit.xyz}}, \code{\link{is.gap}}, \code{\link{gap.inspect}},
\code{\link{orient.pdb}}, \code{\link{pairwise}}, \code{\link{plot.bio3d}}, \code{\link{plot.nma}}, \code{\link{plot.blast}}, \code{\link{biounit}}, etc.).
}
\note{
The latest version, package vignettes and documentation with worked example
outputs can be obtained from the bio3d website:\cr
\url{http://thegrantlab.org/bio3d/}.\cr
\url{http://thegrantlab.org/bio3d/html/}.\cr
\url{http://bitbucket.org/Grantlab/bio3d}.
}
\author{
Barry Grant <bjgrant@umich.edu>
Xin-Qiu Yao <xinqyao@umich.edu>
Lars Skjaerven <larsss@gmail.com>
Julien Ide <julien.ide.fr@gmail.com>
}
\references{
Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
Skjaerven, L. et al. (2014) \emph{BMC Bioinformatics} \bold{15}, 399.
}
\keyword{ documentation }
\examples{
help(package="bio3d")     # list the functions within the package
#lbio3d()                  # list bio3d function names only

## Or visit:
##   http://thegrantlab.org/bio3d/html/

## See the individual functions for further documentation and examples, e.g.
#help(read.pdb)

## Or online:
##    http://thegrantlab.org/bio3d/html/read.pdb.html

\dontrun{
##-- See the list of Bio3D demos
demo(package="bio3d")

## Try some out, e.g:
demo(pdb) # PDB Reading, Manipulation, Searching and Alignment
demo(pca) # Principal Component Analysis
demo(md)  # Molecular Dynamics Trajectory Analysis
demo(nma) # Normal Mode Analysis

## See package vignettes and tutorals online:
##   http://thegrantlab.org/bio3d/tutorials
}
}
