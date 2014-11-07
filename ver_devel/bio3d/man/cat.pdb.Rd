\name{cat.pdb}
\alias{cat.pdb}
\title{ Concatenate two PDB Objects. }
\description{
  Produce a new concatenated PDB object from two smaller PDB objects.
}
\usage{
cat.pdb(pdb1, pdb2, renumber=TRUE, rechain=TRUE)
}
\arguments{
  \item{pdb1}{ a PDB structure object obtained from
    \code{\link{read.pdb}}. }
  \item{pdb2}{ a PDB structure object obtained from
    \code{\link{read.pdb}}. }
  \item{renumber}{ logical, if \sQuote{TRUE} residues will be
    renumbered. }
  \item{rechain}{ logical, if \sQuote{TRUE} molecules will be assigned
    new chain identifiers. }
}
\details{
  This is a basic utility function for creating a concatenated PDB
  object based on two smaller PDB objects. 
}
\value{
  Returns an object of class \code{"pdb"}. See \code{\link{read.pdb}} for
  further details. 
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Lars Skjaerven }
\seealso{ \code{\link{read.pdb}}, \code{\link{atom.select}},
  \code{\link{write.pdb}}, \code{\link{trim.pdb}} }
\examples{
\dontrun{
## Read a PDB file from the RCSB online database
pdb1 <- read.pdb("1etl")
pdb2 <- read.pdb("1hel")

## Concat
new.pdb <- cat.pdb(pdb1, pdb2)

## Write to file
write.pdb(new.pdb, file="concat.pdb")
}
}
\keyword{ utilities }
