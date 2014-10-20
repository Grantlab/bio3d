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
pdb <- read.pdb( system.file("examples/1hel.pdb", package="bio3d") )

# Print structure summary
atom.select(pdb)

# Select all C-alpha atoms with residues numbers between 43 and 54
ca.inds <- atom.select(pdb, "calpha", resno=43:54)
is.select(ca.inds)

}
\keyword{ classes }
