\name{rmsd}
\alias{rmsd}
\title{ Root Mean Square Deviation }
\description{
  Calculate the RMSD between coordinate sets.
}
\usage{
rmsd(a, b=NULL, a.inds=NULL, b.inds=NULL, fit=FALSE, ncore=1, nseg.scale=1)
}
\arguments{
  \item{a}{ a numeric vector containing the reference coordinate set for
    comparison with the coordinates in \code{b}.  Alternatively, if
    \code{b=NULL} then \code{a} can be a matrix or list object
    containing multiple coordinates for pairwise comparison. }
  \item{b}{ a numeric vector, matrix or list object with an \code{xyz}
    component, containing one or more coordinate sets to be compared with
    \code{a}. }
  \item{a.inds }{ a vector of indices that selects the elements of
    \code{a} upon which the calculation should be based. }
  \item{b.inds }{ a vector of indices that selects the elements of
    \code{b} upon which the calculation should be based. }
  \item{fit }{logical, if TRUE coordinate superposition is performed
    prior to RMSD calculation. }  
  \item{ncore }{ number of CPU cores used to do the calculation.
    \code{ncore>1} requires package \sQuote{parallel} installed. } 
  \item{nseg.scale }{ split input data into specified number of segments
    prior to running multiple core calculation. See \code{\link{fit.xyz}}. }
}
\details{
  RMSD is a standard measure of structural distance between
  coordinate sets.

  Structure \code{a[a.inds]} and \code{b[b.inds]} should have the
  same length.

  A least-squares fit is performed prior to RMSD calculation by setting
  \code{fit=TRUE}.  See the function \code{fit.xyz} for more
  details of the fitting process. 
}
\value{
  Returns a numeric vector of RMSD value(s).
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Barry Grant }
\seealso{ \code{\link{fit.xyz}}, \code{\link{rot.lsq}},
  \code{\link{read.pdb}}, \code{\link{read.fasta.pdb}} }
\examples{
# -- Calculate RMSD between two or more structures
aln <- read.fasta(system.file("examples/kif1a.fa",package="bio3d"))
pdbs <- read.fasta.pdb(aln)

# Indices 
gaps <- unique( which(is.na(pdbs$xyz),arr.ind=TRUE)[,2] )
inds <- c(1:ncol(pdbs$xyz))[-gaps]

# RMSD between structure 1 and structure 2
rmsd(a=pdbs$xyz[1,], b=pdbs$xyz[2,], a.inds=inds, b.inds=inds)
rmsd(a=pdbs$xyz[1,], b=pdbs$xyz[2,], a.inds=inds, b.inds=inds, fit=TRUE)

\dontrun{
# RMSD between structure 1 and structures 2 and 3
rmsd(a=pdbs$xyz[1,], b=pdbs$xyz[2:3,], a.inds=inds, b.inds=inds)

# RMSD between structure 1 and all structures in alignment
rmsd(a=pdbs$xyz[1,], b=pdbs, a.inds=inds, b.inds=inds, fit=TRUE)

# pairwise RMSD
rmsd(pdbs$xyz[,inds])
}
}
\keyword{ utilities }
