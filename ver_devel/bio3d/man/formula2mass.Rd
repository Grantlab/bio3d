\name{formula2mass}
\alias{formula2mass}
\title{Chemical Formula to Mass Converter}
\description{Compute the molar mass associated to a chemical formula.}
\usage{
formula2mass(form, sum.mass = TRUE)
}
\arguments{
  \item{form}{a character string containing a chemical formula on the form: 'C3 H5 N O1'.}
  \item{sum.mass}{logical, should the mass of each element be summed.}
}
\details{
Compute the molar mass (in g.mol-1) associated to a chemical formula.
}
\value{
Return a single element numeric vector containing the mass corresponding to a given chemical formula.
}

\author{Lars Skjaerven}

\seealso{
\code{\link{atom2ele}}, \code{\link{atom2mass}}
}
\examples{
#formula2mass("C5 H6 N O3")
}
\keyword{ utilities }
