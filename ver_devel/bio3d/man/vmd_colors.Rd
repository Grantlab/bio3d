\name{vmd_colors}
\alias{vmd_colors}
\title{ VMD Color Palette }
\description{
  This function creates a character vector of the colors used by the VMD
  molecular graphics program.
}
\usage{
vmd_colors(n=33, picker=FALSE, ...)
}
\arguments{
  \item{n}{ The number of desired colors chosen in sequence
    from the VMD color palette (>=1) }
  \item{picker}{ Logical, if TRUE a color wheel plot will be produced to
    aid with color choice. }
  \item{\dots}{ Extra arguments passed to the \code{rgb} function, including
    alpha transparency. }
}
\details{
  The function uses the underlying 33 RGB color codes from VMD, See
  \url{http://www.ks.uiuc.edu/Research/vmd/}. Note that colors will be
  recycled if \dQuote{n} > 33 with a warning issued. When \sQuote{picker} is
  set to \dQuote{TRUE} a color wheel of the requested colors will be
  plotted to the currently active device.
}
\value{
  Returns a character vector with color names.
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.

  \url{http://www.ks.uiuc.edu/Research/vmd/}
}
\author{
  Barry Grant
}
\seealso{
  \code{\link{bwr.colors}}
}
\examples{
## Generate a vector of 10 colors
clrs <- vmd_colors(10)
vmd_colors(4, picker=TRUE)
}
\keyword{utility}
