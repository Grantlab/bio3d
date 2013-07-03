\name{pdbfit}
\alias{pdbfit}
\title{ PDB File Coordinate Superposition }
\description{ Protein Databank Bank file coordinate superposition with
  the Kabsch algorithm. }
\usage{
pdbfit(pdbs, inds = NULL, outpath = NULL, ...)
}
\arguments{
  \item{pdbs}{ a list of class \code{"3dalign"} containing PDB file
    data, as obtained from \code{read.fasta.pdb} or
    \code{pdbaln}. }
  \item{inds}{ a vector of indices that selects the coordinate
    positions, in terms of x, y and z elements, upon which fitting
    should be based. This defults to all equivalent non-gap positions. }
  \item{outpath}{ character string specifing the output directory for
    optional coordinate file output. Note that full files (i.e. all atom
    files) are written, seebelow. }
  \item{\dots}{ extra arguments passed to \code{fit.xyz} function. }
}
\details{
  The function \code{pdbfit} is a wrapper for the function
  \code{fit.xyz}, wherein full details of the superposition procedure
  are documented.

  Input to this function should be a list object obtained with the
  function \code{\link{read.fasta.pdb}} or  \code{\link{pdbaln}}.  See
  the examples below.

  The reference frame for supperposition (i.e. the fixed structure to
  which others are superposed) is the first entry in the input
  \code{"pdbs"} object. For finer control use \code{\link{fit.xyz}}.
}
\value{
  Returns moved coordinates.
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.

  Kabsch \emph{Acta Cryst} (1978) \bold{A34}, 827--828.
}
\author{ Barry Grant }
\seealso{ \code{\link{pdbaln}}, \code{\link{read.fasta.pdb}},
  \code{\link{fit.xyz}}, \code{\link{rmsd}}, \code{\link{read.pdb}} }
\examples{
\dontrun{
#files <- get.pdb(c("4q21","5p21"), URLonly=TRUE)
files <- get.pdb(c("4q21","5p21"), path=tempdir())
pdbs <- pdbaln(files)
xyz <- pdbfit(pdbs)

# Superpose again this time outputing PDBs
xyz <- pdbaln( pdbs, outpath="fitted/" )
}
}
\keyword{ utilities }
