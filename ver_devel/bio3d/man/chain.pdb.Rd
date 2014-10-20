\name{chain.pdb}
\alias{chain.pdb}
\title{ Find Possible PDB Chain Breaks }
\description{
  Find possible chain breaks based on connective Calpha atom separation.
}
\usage{
chain.pdb(pdb, ca.dist = 4, blank = "X", rtn.vec = TRUE)
}
\arguments{
  \item{pdb}{ a PDB structure object obtained from
    \code{\link{read.pdb}}. }
  \item{ca.dist}{ the maximum distance that separates Calpha atoms
    considered to be in the same chain. }
  \item{blank}{ a character to assign non-protein atoms. }
  \item{rtn.vec}{ logical, if TRUE then the one-letter chain vector
    consisting of the 26 upper-case letters of the Roman alphabet is
    returned. }
}
\details{
  This is a basic function for finding possible chain breaks in PDB
  structure files, i.e. connective Calpha atoms that are further than
  \code{ca.dist} apart.
}
\value{
  Prints basic chain information and if \code{rtn.vec} is TRUE returns a
  character vector of chain ids consisting of the 26 upper-case letters of the
  Roman alphabet plus possible \code{blank} entries for non-protein atoms.
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}  
\author{ Barry Grant }
\seealso{ \code{\link{read.pdb}}, \code{\link{atom.select}}, \code{\link{trim.pdb}},
  \code{\link{write.pdb}} }
\examples{
\donttest{
full.pdb <- read.pdb( get.pdb("5p21", URLonly=TRUE) )
inds <- atom.select(full.pdb, resno=c(10:20,30:33))
cut.pdb <- trim.pdb(full.pdb, inds)
chain.pdb(cut.pdb)
}
}
\keyword{ utilities }
