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
nma(pdb, inds = NULL, ff = 'calpha', pfc.fun = NULL, fc.weights = NULL,
    mass = TRUE, temp = 300.0, keep = NULL, cutoff = 15, gamma = 1,
    ncore = 1 )
build.hessian(xyz, pfc.fun, mass.weights = NULL, fc.weights = NULL,
    ncore = 1 ) 
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
  \item{pfc.fun}{ customized pair force constant (\sQuote{pfc})
    function. The provided function should take a vector of distances as
    an argument to return a vector of force constants. If provided,
    'pfc.fun' will override argument \code{ff}. }
  \item{fc.weights}{ a numeric matrix of size NxN (where N is the number
    of calpha atoms) containg scaling factors for the pariwise force
    constants. }
  \item{mass}{ logical, if TRUE the Hessian will be mass-weighted. }
  \item{temp}{ numerical, temperature for which the amplitudes for
    scaling the atomic displacement vectors are calculated. Set
    \sQuote{temp=NULL} to avoid scaling. }
  \item{keep}{ numerical, final number of modes to be stored. Note that
    all subsequent analyses are limited to this subset of modes. }
  \item{ncore}{ number of CPU cores used to build the Hessian matrix.
    Requires multicore package installed. }
  \item{cutoff}{ numerical, cutoff for pair-wise interactions (applies to
    \sQuote{anm}-FF only). }
  \item{gamma}{ numerical, global scaling factor  (applies to
    \sQuote{anm}-FF only). }
  \item{xyz}{ a numeric vector of Cartesian coordinates. }
  \item{mass.weights}{ a numeric vector of atom weights. }
  \item{x}{ an \code{nma} object obtained from \code{\link{nma}}. }
  \item{nmodes}{ numeric, number of modes to be printed. }
  \item{...}{ additional arguments to \sQuote{print}. }
}
\details{
  This function calculates the normal modes of a C-alpha model of a
  protein structure. Two force fields are implemented - both employing
  the elastic network model (ENM).

  The \sQuote{calpha} force field - originally developed by Konrad
  Hinsen - is the recommended one for most applications. It employs a
  spring force constant differentiating between nearest-neighbour pairs
  along the backbone and all other pairs. The force constant function
  was parameterized by fitting to a local minimum of a crambin model
  using the AMBER94 force field. 

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
  \code{mass=FALSE}, and \code{temp=NULL}.
}
\value{
  Returns an object of class \sQuote{nma} with the following components:
  \item{modes}{ numeric matrix  with columns containing the normal mode
    vectors. Mode vectors are converted to unweighted Cartesian
    coordinates  when \code{mass=TRUE}. Note that the 6 first trivial
    eigenvectos appear in columns one to six. }
  \item{frequencies}{ numeric vector containing the vibrational
    frequencies corresponding to each mode (for \code{mass=TRUE}). }
  \item{force.constants}{ numeric vector containing the force constants
    corresponding to each mode (for \code{mass=FALSE)}). }
  \item{fluctuations}{ numeric vector of atomic fluctuations. }
  \item{U}{ numeric vector containing the raw eigenvectors. Equals to
    the \code{modes} component when \code{mass=FALSE} and
    \code{temp=NULL}. }
  \item{L}{ numeric vector containing the raw eigenvalues. }
  \item{xyz}{ numeric vector of the cartesian coordinates in which the
    calculation was performed. }
  \item{mass}{ numeric vector containing the residue masses used for the 
    mass-weighting. }
  \item{temp}{ numerical, temperature for which the amplitudes for
    scaling the atomic displacement vectors are calculated. }
  \item{triv.modes}{ number of trivial modes. }
  \item{natoms}{ number of C-alpha atoms. }
  \item{call}{ the matched call. }
}
\note{
  The current version provides an efficent implementation of NMA with
  execution time comparable to similar software (when the entire Hessian
  is diagonalized).

  The main (speed related) bottleneck is currently the diagonalization
  of the Hessian matrix which is performed with the core R function
  \sQuote{eigen}. For computing a few (5-20) approximate modes the user
  can consult package \sQuote{irlba}.

  NMA is memory extensive and users should be cautions when running
  larger proteins (>3000 residues). Use \sQuote{keep} to reduce the
  amount of memory needed to store the final \sQuote{nma} object (the
  full 3Nx3N Hessian matrix still needs to be allocated).
}
\references{
  Hinsen, K. et al. (2000) \emph{Chemical Physics} \bold{261}, 25--37.
  Atilgan, A.R. et al. (2001) \emph{Biophysical Journal} \bold{80}, 505--515.
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Lars Skjaerven }
\seealso{ \code{\link{fluct.nma}}, \code{\link{mktrj.nma}},
  \code{\link{dccm.nma}}, \code{\link{overlap}}, \code{\link{rmsip}}  }
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
modes <- nma(pdb, ff="anm", mass=FALSE, temp=NULL)

## Customize energy function
## Note - r is a vector of distances !
"my.ff" <- function(r) {
    return(1/r)
}

## Modes with a new funny energy function
modes <- nma(pdb, pfc.fun=my.ff)

## A more manual approach
"anm.ff" <- function(r) {
  ifelse( r>15, 0, 1 )
}

sele <- atom.select(pdb, "//A////CA/")
xyz <- pdb$xyz[sele$xyz]
hessian <- build.hessian(xyz, anm.ff)
modes <- eigen(hessian)

}
\keyword{ analysis }
