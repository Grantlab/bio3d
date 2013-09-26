\name{dccm.pymol}
\alias{dccm.pymol}
\title{ Visualization of Dynamic Cross-Correlation }
\description{
  Visualize a cross-correlation matrix in PyMOL. 
}
\usage{
dccm.pymol(dccm, xyz, step=0.2, omit=0.2)
}
\arguments{
  \item{dccm}{ an object of class \code{dccm} as obtained from
    function \code{dccm} or \code{dccm.nma}. }
  \item{xyz}{ numerical vector, Cartesian coordinates. }
  \item{step}{ binning interval of cross-correlation coefficents. }
  \item{omit}{ correlation coefficents with values (0-omit, 0+omit) will
    be omitted from visualization. }
}
\details{
  This function generates a PyMOL (python) script that will draw colored
  lines between (anti)correlated residues.
}
\value{
  Called for its action.
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Lars Skjaerven }
\seealso{ \code{\link{nma}}, \code{\link{dccm}} }
\examples{
\dontrun{
## Fetch stucture
pdb <- read.pdb("1hel")

## Calculate normal modes
modes <- nma(pdb)

## Calculate correlation matrix
cm <- dccm.nma(modes)

dccm.pymol(cm, modes$xyz)
}
}
\keyword{ utilities }