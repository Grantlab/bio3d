\name{bwr.colors}
\alias{bwr.colors}
\alias{mono.colors}
\title{ Color Palettes }
\description{
  Create a vector of \sQuote{n} \dQuote{contiguous} colors forming
  either a Blue-White-Red or a White-Gray-Black color palette.
}
\usage{
bwr.colors(n)
mono.colors(n)
}
\arguments{
  \item{n}{ the number of colors in the palette (>=1). }
}
\details{
  The function \code{bwr.colors} returns a vector of \code{n} color
  names that range from blue through white to red.
  
  The function \code{mono.colors} returns color names ranging from
  white to black.  Note: the first element of the returned vector
  will be NA. 
}
\value{
  Returns a character vector, \code{cv}, of color names.  This can be
  used either to create a user-defined color palette for subsequent
  graphics with \code{palette(cv)}, or as a \code{col=} specification in
  graphics functions and \code{par}.
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.

  The \code{bwr.colors} function is derived from the \code{gplots}
  package function \code{colorpanel} by Gregory R. Warnes.
}
\author{ Barry Grant }
\seealso{ \code{\link{vmd.colors}}, \code{\link{cm.colors}},
  \code{\link{colors}}, \code{\link{palette}}, \code{\link{hsv}},
  \code{\link{rgb}}, \code{\link{gray}}, \code{\link{col2rgb}}
}
\examples{
\donttest{

# Color a distance matrix
pdb <- read.pdb( "1bg2" )
d <- dm(pdb,"calpha")

plot(d, color.palette=bwr.colors)

plot(d,
     resnum.1 = pdb$atom[pdb$calpha,"resno"],
     color.palette = mono.colors,
     xlab="Residue Number", ylab="Residue Number")

}
}
\keyword{ utilities }

