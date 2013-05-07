\name{rgyr}
\alias{rgyr}
\title{ Radius of Gyration}
\description{
  Calculate the radius of gyration of coordinate sets.
}
\usage{
rgyr(xyz, mass=NULL, ncore=1, nseg.scale=1)
}
\arguments{
  \item{xyz}{ a numeric vector, matrix or list object with an \code{xyz}
    component, containing one or more coordinate sets.}
  \item{mass}{ a numeric vector of atomic masses (unit a.m.u.),
    or a PDB object with masses stored in the "B-factor" column. 
    If \code{mass==NULL}, all atoms are assumed carbon.}
  \item{ncore }{ number of CPU cores used to do the calculation.
    \code{ncore>1} requires multicore package installed. } 
  \item{nseg.scale }{ split input data into specified number of segments
    prior to running multiple core calculation. See \code{\link{fit.xyz}}. }
}
\details{
  Radius of gyration is a standard measure of overall structural change of
  macromolecules.
}
\value{
  Returns a numeric vector of radius of gyration.
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Xin-Qiu Yao & Pete Kekenes-Huskey}
\seealso{ \code{\link{fit.xyz}}, \code{\link{rmsd}},
  \code{\link{read.pdb}}, \code{\link{read.fasta.pdb}} }
\examples{
# -- Calculate Rog of single structure
pdb <- read.pdb("1bg2")
mass <- rep(12, length(pdb$xyz)/3)
mass[substr(pdb$atom[,"elety"], 1, 1) == "N"] <- 14
mass[substr(pdb$atom[,"elety"], 1, 1) == "H"] <- 1
mass[substr(pdb$atom[,"elety"], 1, 1) == "O"] <- 16 
mass[substr(pdb$atom[,"elety"], 1, 1) == "S"] <- 32

rgyr(pdb, mass)

\dontrun{
# -- Calculate Rog of a trajectory
xyz <- read.dcd(system.file("examples/hivp.dcd", package="bio3d"))
rg <- rgyr(xyz)
rg[1:10]

}
}
\keyword{ utilities }
