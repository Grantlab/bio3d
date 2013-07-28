\name{nma}
\alias{nma}
\alias{build.hessian}
\alias{print.nma}
\title{ Normal Modes Analysis }
\description{
  Perform elastic network model (ENM) C-alpha normal modes calculation
  of a protein structure.
}
\usage{
nma(pdb, inds = NULL, ff = 'calpha', pfc.fun = NULL,
    normalize = TRUE, mass = TRUE, temp = 300.0, compiler = TRUE,
    cutoff = 15, gamma = 1)
build.hessian(xyz, pfc.fun, normalize=TRUE)
\method{print}{nma}(x, nmodes=6, \dots)
}
\arguments{
  \item{pdb}{ an object of class \code{pdb} as obtained from
    function \code{read.pdb}. }
  \item{inds}{ atom and xyz coordinate indices obtained from \code{atom.select}
    that selects the elements of \code{pdb} upon which the calculation
    should be based.}
  \item{ff}{ character string specifying the force field to use:
    \sQuote{calpha} or \sQuote{anm}. }
  \item{pfc.fun}{ customized pair force constant function. }
  \item{normalize}{ logical, if TRUE distance vectors will be normalized. }
  \item{mass}{ logical, if TRUE the calculation will be mass
    weighted. }
  \item{temp}{ numerical, temperature at which the atomic vectors are
    calculated. }
  \item{compiler}{ logical, if TRUE code for building the Hessian matrix
    will be byte compiled for efficiency. }
  \item{cutoff}{ numerical, cutoff for pair-wise interactions (applies to
    \sQuote{anm}-FF only). }
  \item{gamma}{ numerical, global scaling factor  (applies to
    \sQuote{anm}-FF only). }
  \item{xyz}{ a numeric vector of Cartesian coordinates. }
  \item{x}{ an \code{nma} object obtained from \code{\link{nma}}. }
  \item{nmodes}{ numeric, number of modes to be printed. }
  \item{...}{ additional arguments to \sQuote{print}. }
}
\details{
  This function calculates the normal modes of a C-alpha model of a
  protein structure. Two force fields are implemented - both employing
  the elastic network model (ENM).

  The \sQuote{Calpha} force field - originally developed by Konrad
  Hinsen - is the recommended one for most applications. It employs a
  spring force constant differentiating between nearest-neighbour pairs
  along the backbone and all other pairs. The force constant is
  obtained by fitting to a local minimum of the AMBER94 force field.

  The implementation of the \sQuote{ANM} (Anisotropic Network Model)
  force field originates from the lab of Ivet Bahar. It uses a
  simplified spring force constant based on the pair-wise
  distance.

  By default \code{nma()} will diagonalize the mass-weighted Hessian
  matrix. The resulting mode vectors are moreover scaled by the thermal
  fluctuation amplitudes.

  The implementation under default arguments reproduces the calculation
  of normal modes (VibrationalModes) in the Molecular Modeling Toolkit
  (MMTK) package. To reproduce ANM modes set \code{ff='anm'},
  \code{mass=FALSE}, \code{temp=NULL}, and \code{normalize=FALSE}.
}
\value{
  Returns an object of class \sQuote{nma} with the following components:
  \item{modes}{ numeric matrix containing the normal mode vectors
    (converted to unweighted Cartesian coordinates). }
  \item{frequencies}{ numeric vector containing the vibrational
    frequencies corresponding to each mode (for \code{mass=TRUE}). }
  \item{force.constants}{ numeric vector containing the force constants
    corresponding to each mode (for \code{mass=FALSE)}). }
  \item{fluctuations}{ numeric vector of atomic fluctuations . }
  \item{U}{ numeric vector containing the raw eigenvectors. }
  \item{L}{ numeric vector containing the raw eigenvalues. }
  \item{xyz}{ numeric vector of the cartesian coordinates in which the
    calculation was performed. }
  \item{mass}{ numeric vector containing the residue masses used for the 
    mass-weighting. }
  \item{temp}{ numerical, temperature at which the atomic vectors are
    calculated. }
  \item{triv.modes}{ number of trivial modes. }
  \item{natoms}{ number of C-alpha atoms. }
  \item{call}{ the matched call. }
}
\references{
  Hinsen, K. et al. (2000) \emph{Chemical Physics} \bold{261}, 25--37.
  Atilgan, A.R. et al. (2001) \emph{Biophysical Journal} \bold{80}, 505--515.
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Lars Skjaerven }
\seealso{ \code{\link{fluct.nma}}, \code{\link{mktrj.nma}},
  \code{\link{dccm.nma}}, \code{\link{overlap}}  }
\examples{
## Fetch stucture
pdb <- read.pdb("1hel")

## Calculate normal modes
modes <- nma(pdb)

## Print modes
print(modes)

## Plot modes
plot(modes)

## Visualize modes
m7 <- mktrj.nma(modes, mode=7, file="mode_7.pdb")

## Calculate 'energetic' modes
modes <- nma(pdb, mass=FALSE)

## Use Anisotropic Network Model
modes <- nma(pdb, ff="anm", mass=FALSE, temp=NULL, normalize=FALSE)

## Customize energy function
"my.ff" <- function(r) {
    return(1 / (r**2))
}

## Modes with a new funny energy function
modes <- nma(pdb, pfc.fun=my.ff)

}
\keyword{ analysis }
