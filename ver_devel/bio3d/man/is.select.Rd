\name{is.select}
\alias{is.select}
\title{Is an Object of Class \sQuote{select}?}
\description{
Checks whether its argument is an object of class \sQuote{select}.
}
\usage{
is.select(x)
}
\arguments{
  \item{x}{an R object to be tested.}
}
\details{
Tests if x is an object of class \sQuote{select}, i.e. if x has a \dQuote{class} attribute equal to \code{select}.
}
\value{
TRUE if x is an object of class \sQuote{select} and FALSE otherwise
}

\author{
Julien Ide
}

\seealso{
\code{\link{atom.select}}
}
\examples{
# Read a PDB file
pdb <- read.pdb( "http://www.rcsb.org/pdb/files/1BG2.pdb" )

# Print structure summary
atom.select(pdb)

# Select all C-alpha atoms with residues numbers between 65 and 143
ca.inds <- atom.select(pdb, resno=65:143, elety="CA")
is.select(ca.inds)

}
\keyword{ classes }
