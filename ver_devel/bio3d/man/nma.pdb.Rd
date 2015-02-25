\name{nma.pdb}
\alias{nma.pdb}
\alias{build.hessian}
\alias{print.nma}
\title{ Normal Mode Analysis }
\description{
  Perform elastic network model (ENM) C-alpha normal modes calculation
  of a protein structure.
}
\usage{
\method{nma}{pdb}(pdb, inds = NULL, ff = 'calpha', pfc.fun = NULL,
                  mass = TRUE, temp = 300.0, keep = NULL, hessian = NULL,
                  outmodes = NULL, \dots )

build.hessian(xyz, pfc.fun, fc.weights = NULL, sequ = NULL, sse = NULL,
    ss.bonds = NULL, \dots)

\method{print}{nma}(x, nmodes=6, \dots)
}
\arguments{
  \item{pdb}{ an object of class \code{pdb} as obtained from
    function \code{\link{read.pdb}}. }
  \item{inds}{ atom and xyz coordinate indices obtained from
    \code{\link{atom.select}} that selects the elements of \code{pdb}
    upon which the calculation should be based. If not provided the
    function will attempt to select the calpha atoms automatically
    (based on function \code{\link{atom.select}}). }
  \item{ff}{ character string specifying the force field to use:
    \sQuote{calpha}, \sQuote{anm}, \sQuote{pfanm}, \sQuote{calphax},
    \sQuote{reach}, or \sQuote{sdenm}. }
  \item{pfc.fun}{ customized pair force constant (\sQuote{pfc})
    function. The provided function should take a vector of distances as
    an argument to return a vector of force constants. If provided,
    'pfc.fun' will override argument \code{ff}. See examples below. }
  \item{mass}{ logical, if TRUE the Hessian will be mass-weighted. }
  \item{temp}{ numerical, temperature for which the amplitudes for
    scaling the atomic displacement vectors are calculated. Set
    \sQuote{temp=NULL} to avoid scaling. }
  \item{keep}{ numerical, final number of modes to be stored. Note that
    all subsequent analyses are limited to this subset of modes. This
    option is useful for very large structures and cases where memory
    may be limiting. }
  \item{hessian}{ hessian matrix as obtained from
    \code{\link{build.hessian}}. For internal purposes and generally not
    intended for public use. }
  \item{outmodes}{ atom indices as obtained from \code{\link{atom.select}}) 
    specifying the atoms to include in the resulting mode object. }
  \item{xyz}{ a numeric vector of Cartesian coordinates. }
  \item{fc.weights}{ a numeric matrix of size NxN (where N is the number
    of calpha atoms) containg scaling factors for the pariwise force
    constants. See examples below. }
  \item{sse}{ secondary structure elements as obtained from
    \code{dssp}. }
  \item{sequ}{ a character vector of the amino acid sequence. }
  \item{ss.bonds}{ a numeric two-column matrix containing the residue
    numbers of the disulfide bridges in the structure. }
  \item{x}{ an \code{nma} object obtained from \code{\link{nma.pdb}}. }
  \item{nmodes}{ numeric, number of modes to be printed. }
  \item{...}{ additional arguments to \code{\link{build.hessian}},
    \code{\link{aa2mass}}, \code{pfc.fun}, and \code{\link{print}}. One
    useful option here for dealing with unconventional residues is
    \sQuote{mass.custom}, see the \code{\link{aa2mass}} function for details. }
}
\details{
  This function calculates the normal modes of a C-alpha model of a
  protein structure. A number of force fields are implemented all of
  whhich employ the elastic network model (ENM).

  The \sQuote{calpha} force field - originally developed by Konrad
  Hinsen - is the recommended one for most applications. It employs a
  spring force constant differentiating between nearest-neighbour pairs
  along the backbone and all other pairs. The force constant function
  was parameterized by fitting to a local minimum of a crambin model
  using the AMBER94 force field.

  See \code{\link{load.enmff}} for details of the different force fields.

  By default \code{\link{nma.pdb}} will diagonalize the mass-weighted Hessian
  matrix. The resulting mode vectors are moreover scaled by the thermal
  fluctuation amplitudes.

  The implementation under default arguments reproduces the calculation
  of normal modes (VibrationalModes) in the Molecular Modeling Toolkit
  (MMTK) package. To reproduce ANM modes set \code{ff='anm'},
  \code{mass=FALSE}, and \code{temp=NULL}.
}
\value{
  Returns an object of class \sQuote{nma} with the following components:
  \item{modes}{ numeric matrix with columns containing the normal mode
    vectors. Mode vectors are converted to unweighted Cartesian
    coordinates  when \code{mass=TRUE}. Note that the 6 first trivial
    eigenvectos appear in columns one to six. }
  \item{frequencies}{ numeric vector containing the vibrational
    frequencies corresponding to each mode (for \code{mass=TRUE}). }
  \item{force.constants}{ numeric vector containing the force constants
    corresponding to each mode (for \code{mass=FALSE)}). }
  \item{fluctuations}{ numeric vector of atomic fluctuations. }
  \item{U}{ numeric matrix with columns containing the raw
    eigenvectors. Equals to the \code{modes} component when
    \code{mass=FALSE} and \code{temp=NULL}. }
  \item{L}{ numeric vector containing the raw eigenvalues. }
  \item{xyz}{ numeric matrix of class \code{xyz} containing the
    Cartesian coordinates in which the calculation was performed. }
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
  \code{\link{eigen}}. For computing a few (5-20) approximate modes the user
  can consult package \sQuote{irlba}.

  NMA is memory extensive and users should be cautions when running
  larger proteins (>3000 residues). Use \sQuote{keep} to reduce the
  amount of memory needed to store the final \sQuote{nma} object (the
  full 3Nx3N Hessian matrix still needs to be allocated).

  We thank Edvin Fuglebakk for valuable discussions on the
  implementation as well as for contributing with testing.
}
\references{
  Skjaerven, L. et al. (2014) \emph{BMC Bioinformatics} \bold{15}, 399.
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
  Hinsen, K. et al. (2000) \emph{Chemical Physics} \bold{261}, 25--37.
}
\author{ Lars Skjaerven }
\seealso{ \code{\link{fluct.nma}}, \code{\link{mktrj.nma}},
  \code{\link{dccm.nma}}, \code{\link{overlap}}, \code{\link{rmsip}},
  \code{\link{load.enmff}}. }
\examples{
## Fetch stucture
pdb <- read.pdb( system.file("examples/1hel.pdb", package="bio3d") )

## Calculate normal modes
modes <- nma(pdb)

## Print modes
print(modes)

## Plot modes
plot(modes)

## Visualize modes
#m7 <- mktrj.nma(modes, mode=7, file="mode_7.pdb")

\dontrun{
## Use Anisotropic Network Model
modes <- nma(pdb, ff="anm", mass=FALSE, temp=NULL, cutoff=15)

## Use SSE information and SS-bonds
sse      <- dssp(pdb, resno=FALSE, full=TRUE)
ss.bonds <- matrix(c(76,94, 64,80, 30,115, 6,127),
                   ncol=2, byrow=TRUE)

modes <- nma(pdb, ff="calphax", sse=sse, ss.bonds=ss.bonds)


## User defined energy function
## Note: Must take a vector of distances
"my.ff" <- function(r) {
   ifelse( r>15, 0, 1 )
}

## Modes with a user defined energy function
modes <- nma(pdb, pfc.fun=my.ff)


## A more manual approach
sele <- atom.select(pdb, chain='A', elety='CA')
xyz <- pdb$xyz[sele$xyz]

hessian <- build.hessian(xyz, my.ff)
modes <- eigen(hessian)

## Dealing with unconventional residues
pdb <- read.pdb("1xj0")

## nma(pdb)
modes <- nma(pdb, mass.custom=list(CSX=121.166))
}
}
\keyword{ analysis }
