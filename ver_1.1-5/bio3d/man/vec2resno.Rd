\name{vec2resno}
\alias{vec2resno}
\title{ Replicate Per-residue Vector Values }
\description{
  Replicate values in one vector based on consecutive entries in a second
  vector. Useful for adding per-residue data to all-atom PDB files.
}
\usage{
vec2resno(vec, resno)
}
\arguments{
  \item{vec}{ a vector of values to be replicated.}
  \item{resno}{ a reference vector or a PDB structure object, obtained
    from \code{\link{read.pdb}}, upon which replication is based. }
}
\details{
  This function can aid in mapping data to PDB structure files. For
  example, residue conservation per position (or any other one value per
  residue data) can be replicated to fit the B-factor field of an all
  atom PDB file which can then be rendered according to this field in a
  molecular viewer.
  
  A basic check is made to ensure that the number of consecutively
  unique entries in the reference vector equals the length of the vector
  to be replicated.
}
\value{
  Returns a vector of replicated values.
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}  
\author{ Barry Grant }
\seealso{ \code{\link{read.pdb}}, \code{\link{atom.select}},
  \code{\link{write.pdb}} }
\examples{
vec2resno(c("a","b"), c(1,1,1,1,2,2))
}
\keyword{ utilities }
