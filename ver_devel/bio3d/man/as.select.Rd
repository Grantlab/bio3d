\name{as.select}
\alias{as.select}
\title{ Convert Atomic Indices to a Select Object }
\description{
  Convert atomic indices to a select object with \sQuote{atom} and
  \sQuote{xyz} components.
}
\usage{
as.select(x, \dots)
}
\arguments{
  \item{x}{ a numeric vector containing atomic indices to be converted
    to a \sQuote{select} object. Alternatively, a logical vector can be
    provided. }
  \item{\dots}{ arguments passed to and from functions. }
}
\details{
  Convert atomic indices to a select object with \sQuote{atom} and
  \sQuote{xyz} components. 
}
\value{
  Returns a list of class \code{"select"} with the following components:
  \item{atom}{ a numeric matrix of atomic indices. }
  \item{xyz }{ a numeric matrix of xyz indices. }
  \item{call }{ the matched call. }
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Lars Skjaerven }
\seealso{
  \code{\link{atom.select}}, \code{\link{read.pdb}}
}
\examples{
as.select(c(1,2,3))
}
\keyword{ utilities }
