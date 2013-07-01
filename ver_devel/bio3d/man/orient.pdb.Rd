\name{orient.pdb}
\alias{orient.pdb}
\title{ Orient a PDB Structure }
\description{
   Center, to the coordinate origin, and orient, by principal axes, the
   coordinates of a given PDB structure or xyz vector.
 }
\usage{
orient.pdb(pdb, atom.subset = NULL, verbose = TRUE)
}
\arguments{
  \item{pdb}{ a pdb data structure obtained from \code{\link{read.pdb}}
    or a vector of \sQuote{xyz} coordinates. }
  \item{atom.subset}{ a subset of atom positions to base orientation on. }
  \item{verbose}{ print dimension details. }
}
\value{
  Returns a numeric vector of re-oriented coordinates.
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Barry Grant }
\note{
  Centering and orientation can be restricted to a \code{atom.subset} of atoms.
}
\seealso{ \code{\link{read.pdb}}, \code{\link{write.pdb}},
\code{\link{fit.xyz}}, \code{\link{rot.lsq}} , \code{\link{atom.select}}}
\examples{
pdb <- read.pdb( "1bg2" )
xyz <- orient.pdb(pdb)
#write.pdb(pdb, xyz = xyz, file = "mov1.pdb")


# Based on C-alphas
inds <- atom.select(pdb, "calpha")
xyz  <- orient.pdb(pdb, atom.subset=inds$atom)
#write.pdb(pdb, xyz = xyz, file = "mov2.pdb")


# Based on a central Beta-strand
inds <- atom.select(pdb, "///224:232///CA/")
xyz  <- orient.pdb(pdb, atom.subset=inds$atom)
#write.pdb(pdb, xyz = xyz, file = "mov3.pdb")
}
\keyword{ manip }
