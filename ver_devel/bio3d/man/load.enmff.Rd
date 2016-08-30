\name{load.enmff}
\alias{load.enmff}
\alias{ff.calpha}
\alias{ff.anm}
\alias{ff.pfanm}
\alias{ff.sdenm}
\alias{ff.reach}
\alias{ff.aaenm}
\alias{ff.aaenm2}
\title{ ENM Force Field Loader }
\description{
  Load force field for elastic network normal mode calculation.
}
\usage{
load.enmff(ff = 'calpha')
ff.calpha(r, rmin=2.9, ...)
ff.anm(r, cutoff=15, gamma=1, ...)
ff.pfanm(r, cutoff=NULL, ...)
ff.sdenm(r, atom.id, pdb, ...)
ff.reach(r, atom.id, ...)
ff.aaenm(r, ...)
ff.aaenm2(r, atom.id, pdb, ...)
}
\arguments{
  \item{ff}{ a character string specifying the force field to use:
    \sQuote{calpha}, \sQuote{anm}, \sQuote{pfanm}, \sQuote{calphax},
    \sQuote{reach}, or  \sQuote{sdenm}. }
  \item{r}{ a numeric vector of c-alpha distances. }
  \item{rmin}{ lowest allowed atom-atom distance for the force constant
    calculation. The default of 2.9A is based on an evaluation of 24
    high-resolution X-ray structures (< 1A). }
  \item{cutoff}{ numerical, cutoff for pair-wise interactions. }
  \item{gamma}{ numerical, global scaling factor. }
  \item{atom.id}{ atomic index. }
  \item{pdb}{ a \code{pdb} object as obtained from function \code{read.pdb}. }
  \item{...}{ additional arguments (for technical reasons). }
}
\details{
  This function provides a collection of elastic network model (ENM)
  force fields for normal modes analysis (NMA) of protein structures. It
  returns a function for calculating the residue-residue spring force
  constants.

  The \sQuote{calpha} force field - originally developed by Konrad
  Hinsen - is the recommended one for most applications. It employs a
  spring force constant differentiating between nearest-neighbour pairs
  along the backbone and all other pairs. The force constant function
  was parameterized by fitting to a local minimum of a crambin model
  using the AMBER94 force field.

  The implementation of the \sQuote{ANM} (Anisotropic Network Model)
  force field originates from the lab of Ivet Bahar. It uses a
  simplified (step function) spring force constant based on the
  pair-wise distance. A variant of this from the Jernigan lab is the
  so-called \sQuote{pfANM} (parameter free ANM) with interactions that
  fall off with the square of the distance.
  
  The \sQuote{calphax} force field is an extension of the original
  \sQuote{calpha} force field by Hinsen. In this implementation we have
  included specific force constants for disulfide bridges, helix 1-4
  interactions, and beta sheet bridges.

  The \sQuote{sdENM} (by Dehouck and Mikhailov) employs residue specific
  spring force constants. It has been parameterized through a
  statistical analysis of a total of 1500 NMR ensembles.

  The \sQuote{REACH} force field (by Moritsugu and Smith) is
  parameterized based on variance-covariance matrices obtained from MD
  simulations. It employs force constants that fall off exponentially
  with distance for non-bonded pairs.

  The all-atom ENM force fields (\sQuote{aaenm} and \sQuote{aaenm2}) was
  obtained by fitting to a local energy minimum of a crambin model
  derived from the AMBER99SB force field (same approach as in Hinsen et
  al 2000). It employs a pair force constant function which falls as
  r^-6. \sQuote{aanma2} employs additonally specific force constants for
  covalent and intra-residue atom pairs. See also \code{\link{aanma}}
  for more details.

  See references for more details on the individual force fields.
}
\note{
  The arguments \sQuote{atom.id} and \sQuote{pdb} are used
  from within function \sQuote{build.hessian} for functions that are not
  simply a function of the pair-wise distance. e.g. the force constants
  in the \sQuote{sdENM} model computes the force constants based on a
  function of the residue types and calpha distance. 
}
\value{
  \sQuote{load.enmff} returns a function for calculating the spring 
  force constants. The \sQuote{ff} functions returns a numeric vector of
  residue-residue spring force constants. 
}
\references{
  Skjaerven, L. et al. (2014) \emph{BMC Bioinformatics} \bold{15}, 399.
  Hinsen, K. et al. (2000) \emph{Chemical Physics} \bold{261}, 25--37.
  Atilgan, A.R. et al. (2001) \emph{Biophysical Journal} \bold{80}, 505--515.
  Dehouck Y. & Mikhailov A.S. (2013) \emph{PLoS Comput Biol} \bold{9}:e1003209.
  Moritsugu K. & Smith J.C. (2008) \emph{Biophysical Journal} \bold{95}, 1639--1648.
  Yang, L. et al. (2009) \emph{PNAS} \bold{104}, 12347-52.
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Lars Skjaerven }
\seealso{ \code{\link{nma}}, \code{\link{build.hessian}} }
\examples{
## Load the c-alpha force field
pfc.fun <- load.enmff('calpha')

## Calculate the pair force constant for a set of C-alpha distances
force.constants <- pfc.fun( seq(4,8, by=0.5) )

## Calculate the complete spring force constant matrix
## Fetch PDB
pdb <- read.pdb( system.file("examples/1hel.pdb", package="bio3d") )

## Fetch only c-alpha coordinates
ca.inds <- atom.select(pdb, 'calpha')
xyz <- pdb$xyz[ca.inds$xyz]

## Calculate distance matrix
dists <- dm.xyz(xyz, mask.lower=FALSE)

## all pair-wise spring force constants
fc.matrix <- apply(dists, 1, pfc.fun)
}
\keyword{ utilities }

 
  
  
