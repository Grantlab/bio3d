\name{pymol}
\alias{pymol}
\alias{pymol.pdbs}
\title{ Biomolecular Visualization with PyMOL }
\description{
  Visualize Bio3D structure objects in PyMOL
}
\usage{
pymol(\dots)

\method{pymol}{pdbs}(pdbs, col=NULL, file="R.pdb")
}
\arguments{
  \item{pdbs}{ aligned C-alpha Cartesian coordinates as obtained with 
    \code{\link{read.fasta.pdb}} or \code{\link{pdbaln}}. }
  \item{col}{ a single element character vector specifying the coloring
    of the structures. Options are: \sQuote{index}, \sQuote{index2},
    \sQuote{gaps}, \sQuote{rmsf}, \sQuote{putty}.

    Special cases: Provide a \sQuote{core} object as obtained by
    \code{core.find} to color on the invariant core. Alternatively,
    provide a vector containing the color code for each structure in the
    \sQuote{pdbs} object. }
  \item{file}{ a single element character vector specifying the file name
    of the PyMOL session file. }
}
\details{
  These functions provides a convenient approach for the visualization of
  Bio3D objects in PyMOL. See examples for more details. 
}
\value{
  Called for its action
}
\references{
  Skjaerven, L. et al. (2014) \emph{BMC Bioinformatics} \bold{15}, 399.
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Lars Skjaerven }
\seealso{
  \code{\link{view}}
}
\examples{
\dontrun{
attach(transducin)

# build a pymol session containing all structures in the PDBs object
pymol(pdbs)

# color by invariant core (
# core <- core.find(pdbs)
pymol(pdbs, col=core)

# color by RMSF
pymol(pdbs, col="rmsf")

# color by clustering
rd <- rmsd(pdbs$xyz)
hc <- hclust(as.dist(rd))
grps <- cutree(hc, k=3)
pymol(pdbs, col=grps)
}
}