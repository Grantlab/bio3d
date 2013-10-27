\name{deformation}
\alias{deformation}
\title{ Deformation Analysis }
\description{
  Calculate deformation energies from Normal Mode Analysis.
}
\usage{
deformation(nma, mode.inds = NULL, pfc.fun = NULL, ncore = NULL)
}
\arguments{
  \item{nma}{ a list object of class \code{"nma"} (obtained with
    \code{\link{nma}}).}
  \item{mode.inds}{ a numeric vector of mode indices in which the
    calculation should be based. }
  \item{pfc.fun}{ customized pair force constant (\sQuote{pfc})
    function. The provided function should take a vector of distances as
    an argument to return a vector of force constants. See \code{nma}
    for examples. }
  \item{ncore }{ number of CPU cores used to do the calculation.
    \code{ncore>1} requires multicore package installed. }
}
\details{
  Deformation analysis provides a measure for the amount of local
  flexibility of the protein structure - i.e. atomic motion relative to
  neighbouring atoms. It differs from \sQuote{fluctuations} (e.g. RMSF
  values) which provide amplitudes of the absolute atomic motion.
    
  Deformation energies are calculated based on the \code{nma} object. By
  default the first 20 non-trivial modes are included in the calculation.

  See examples for more details.
}
\value{
  Returns a list with the following components:
  \item{ei }{ numeric matrix containing the energy contribution (E) from
    each atom (i; row-wise) at each mode index (column-wise). }
  \item{sums }{ deformation energies corresponding to each mode. }
}
\references{
  Hinsen, K. (1998) \emph{Proteins} \bold{33}, 417--429.
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Lars Skjaerven }
\seealso{ \code{\link{nma} } }
\examples{
## Fetch stucture
pdb <- read.pdb("1hel")

## Calculate (vibrational) normal modes
modes <- nma(pdb)

## Calculate deformation energies
def.energies <- deformation(modes)

\dontrun{
## Fluctuations of first non-trivial mode
def.energies <- deformation(modes, mode.inds=seq(7, 16))

write.pdb(pdb=NULL, xyz=modes$xyz,
          b=def.energies$ei[,1])
}
}
\keyword{ analysis }