\name{print.xyz}
\alias{print.xyz}
\title{ Printing XYZ coordinates }
\description{
    Print method for objects of class \sQuote{xyz}.
}
\usage{
\method{print}{xyz}(x, ...)
}
\arguments{
  \item{x}{ a \sQuote{xyz} object indicating 3-D coordinates of 
    biological molecules. }
  \item{\dots}{ additional arguments passed to \sQuote{print}. }
}
\value{
  Called for its effect.
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.}
\author{ Barry Grant }
\seealso{
  \code{\link{is.xyz}}, \code{\link{read.ncdf}},
  \code{\link{read.pdb}}, \code{\link{read.dcd}}, \code{\link{fit.xyz}}
}
\examples{

# Read a PDB file
pdb <- read.pdb("1BG2")
print(pdb$xyz)

}
\keyword{ utilities }

