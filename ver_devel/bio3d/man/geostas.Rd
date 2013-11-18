\name{geostas}
\alias{geostas}
\alias{amsm.xyz}
\title{ Domain Identifier: Geometrically Stable Domains }
\description{
  Identifies geometrically stable domains 
}
\usage{
geostas(xyz, amsm = NULL, k = 3, \dots)

amsm.xyz(xyz, ncore = NULL)
}
\arguments{
  \item{xyz}{ numeric matrix of xyz coordinates as obtained e.g. by
    \code{\link{read.ncdf}}, \code{\link{read.dcd}}, or
    \code{\link{mktrj}}. }
  \item{amsm}{ a numeric matrix as obtained by \code{\link{amsm.xyz}}. }
  \item{ncore }{ number of CPU cores used to do the calculation.
    \code{ncore>1} requires multicore package installed. }
  \item{k }{ an integer scalar or vector with the desired number of
    groups. }
  \item{...}{ additional arguments to \code{\link{amsm.xyz}}. }
}
\details{
  
  See examples for more details.
}
\value{
  Returns a list object with the following components:
  \item{grps }{ a numeric vector containing the domain assignment per
    residue. }
  \item{amsm }{ a numeric matrix of atomic movement similarity (AMSM). }
}
\references{
  Romanowska, J. et al. (2012) \emph{JCTC} \bold{8}, 2588--2599.
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Julia Romanowska and Lars Skjaerven }
\seealso{ \code{\link{read.ncdf}}, \code{\link{read.pdb}},
  \code{\link{read.dcd}}, \code{\link{mktrj}}. }
\examples{
## Fetch stucture
pdb <- read.pdb("1hel")

## Calculate (vibrational) normal modes
modes <- nma(pdb)

## Mode trajectories
trj <- rbind(mktrj(modes, mode=7),
             mktrj(modes, mode=8),
             mktrj(modes, mode=9))

## Calculate stable domains
gs <- geostas(trj)

## Make a new PDB file with chain ID corresponding to domains
## write.pdb(pdb=NULL, xyz=modes$xyz, chain=gs$grps)

## gs <- geostas(trj, amsm=gs$amsm, k=5)

\dontrun{
dcdfile <- system.file("examples/hivp.dcd", package = "bio3d")
pdbfile <- system.file("examples/hivp.pdb", package = "bio3d")

dcd <- read.dcd(dcdfile)
pdb <- read.pdb(pdbfile)

gs <- geostas(dcd)
write.pdb(pdb, chain=gs$grps)
}
}
\keyword{ analysis }