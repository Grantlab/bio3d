\name{view.dccm}
\alias{view.dccm}
\title{ Visualization of Dynamic Cross-Correlation }
\description{
  Structural visualization of a cross-correlation matrix. 
}
\usage{
view.dccm(dccm, pdb, step=0.2, omit=0.2, radius=0.15, type="pymol", outprefix="corr",
          launch=FALSE, exefile="pymol")
}
\arguments{
  \item{dccm}{ an object of class \code{dccm} as obtained from
    function \code{dccm} or \code{dccm.nma}. }
  \item{pdb}{ an object of class \code{pdb} as obtained from
    function \code{read.pdb} or a numerical vector of Cartesian
    coordinates. } 
  \item{step}{ binning interval of cross-correlation coefficents. }
  \item{omit}{ correlation coefficents with values (0-omit, 0+omit) will
    be omitted from visualization. }
  \item{radius}{ radius of visualized correlations. }
  \item{type}{ character string specifying the type of visualization:
    \sQuote{pymol} or \sQuote{pdb}. }
  \item{outprefix}{ character string specifying the file prefix. If
    \code{NULL} the temp directory will be used. }
  \item{launch}{ logical, if TRUE PyMol will be launched. }
  \item{exefile}{ file path to the \sQuote{PYMOL} program on your system
    (i.e. how is \sQuote{PYMOL} invoked). }
}
\details{
  This function generates a PyMOL (python) script that will draw colored
  lines between (anti)correlated residues. The PyMOL script file is
  stored in the working directory with filename \dQuote{corr.py}, with
  coordinates in PDB format with filename \dQuote{corr.inpcrd.pdb}.
  PyMOL will only be launched when using argument
  \sQuote{launch=TRUE}. Alternatively a PDB file with CONECT records
  will be generated (when argument \code{type='pdb'}). 
   
  For the PyMOL version, PyMOL CGO objects are generated - each object
  representing a range of correlation values (corresponding to the
  actual correlation values as found in the correlation
  matrix). E.g. the PyMOL object with name \dQuote{cor_-1_-08} would
  display all pairs of correlations with values between -1 and -0.8.
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

view.dccm(cm, modes$xyz)
}
}
\keyword{ utilities }
