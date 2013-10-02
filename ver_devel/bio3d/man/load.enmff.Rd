\name{load.enmff}
\alias{load.enmff}
\alias{ff.calpha}
\alias{ff.calphax}
\alias{ff.anm}
\alias{ff.pfanm}
\alias{ff.sdenm}
\alias{ff.reach}
\title{ ENM Force Field Loader }
\description{
  Load force field for elastic network normal mode calculation.
}
\usage{
load.enmff(ff = 'calpha')
ff.calpha(r, ...)
ff.anm(r, cutoff=15, gamma=1, ...)
ff.pfanm(r, cutoff=NULL, ...)
ff.calphax(r, atom.id, ssdat=NULL, verbose=FALSE, ...)
ff.sdenm(r, atom.id, ssdat=NULL, ...)
ff.reach(r, atom.id, ssdat=NULL, ...)
}
\arguments{
  \item{ff}{ a character string specifying the force field to use:
    \sQuote{calpha}, \sQuote{anm}, \sQuote{pfanm}, \sQuote{calphax},
    \sQuote{reach}, or  \sQuote{sdenm}. }
  \item{r}{ a numeric vector of c-alpha distances. }
  \item{cutoff}{ numerical, cutoff for pair-wise interactions. }
  \item{gamma}{ numerical, global scaling factor. }
  \item{atom.id}{ atomic index. }
  \item{ssdat}{ sequence and structure data. }
  \item{verbose}{ logical, if TRUE interaction details are printed. }
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
  pair-wise distance. A variant of this is the so-called \sQuote{pfANM}
  (parameter free ANM) with interactions that fall off with the square
  of the distance. 
  
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

  See references for more details on the individual force fields.
}
\note{
  The arguments \sQuote{atom.id} and \sQuote{ssdat} are used
  from within function \sQuote{build.hessian} for functions that are not
  simply a function of the pair-wise distance. e.g. the force constants
  in the \sQuote{calphax} model is a function of c-alpha distances, SSE 
  information and SS-bonds, while the \sQuote{sdENM} force field
  computes the force constants based on a function of the residue types
  and calpha distance.
}
\value{
  \sQuote{load.enmff} returns a function for calculating the spring 
  force constants. The \sQuote{ff} functions returns a numeric vector of
  residue-residue spring force constants. 
}
\references{
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
pdb <- read.pdb("1hel")

## Fetch only c-alpha coordinates
ca.inds <- atom.select(pdb, 'calpha')
xyz <- pdb$xyz[ca.inds$xyz]

## Calculate distance matrix
dists <- dm.xyz(xyz, mask.lower=FALSE)

## all pair-wise spring force constants
fc.matrix <- apply(dists, 1, pfc.fun)
}
\keyword{ utilities }

 
  
  