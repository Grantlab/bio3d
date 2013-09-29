\name{view.modes}
\alias{view.modes}
\title{ Vector Field Visualization of Modes }
\description{
  Structural visualization of mode vectors obtained from PCA or NMA.
}
\usage{
view.modes(modes, mode=NULL, outprefix="mode_vecs",
          scale=5, dual=FALSE, launch=FALSE)
}
\arguments{
  \item{modes}{ an object of class \code{nma} or \code{pca} as obtained
    from functions \code{nma} or \code{pca.xyz}. }
  \item{mode}{ the mode number for which the vector field should be
    made. } 
  \item{outprefix}{ character string specifying the file prefix. If
    \code{NULL} the temp directory will be used. }
  \item{scale}{ global scaling factor. }
  \item{dual}{ logical, if TRUE mode vectors are also drawn in both
    direction. }
  \item{launch}{ logical, if TRUE PyMol will be launched. }
}
\details{
  This function generates a PyMOL (python) script for drawing
  mode vectors on a PDB structure.
}
\value{
  Called for its action.
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Lars Skjaerven }
\seealso{ \code{\link{nma}}, \code{\link{pca.xyz}} }
\examples{
\dontrun{
## Fetch stucture
pdb <- read.pdb("1hel")

## Calculate normal modes
modes <- nma(pdb)

view.modes(modes, mode=7)
}
}
\keyword{ utilities }
