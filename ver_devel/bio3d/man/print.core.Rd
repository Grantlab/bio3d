\name{print.core}
\alias{print.core}
\title{ Printing Core Positions and Returning Indices }
\description{
  Print method for core.find objects.
}
\usage{
\method{print}{core}(x, vol = NULL, ...)
}
\arguments{
  \item{x}{ a list object obtained with the function
    \code{\link{core.find}}. }
  \item{vol}{ the maximal cumulative volume value at which core positions
    are detailed. }
  \item{...}{ additional arguments to \sQuote{print}. }
}
\value{
  Returns a three component list of indices:
  \item{atom}{atom indices of core positions}
  \item{xyz}{xyz indices of core positions}
  \item{resno}{residue numbers of core positions}
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Barry Grant }
\note{
  The produced \code{\link{plot.core}} function can be useful for deciding on the
  core/non-core boundary.
}
\seealso{ \code{\link{core.find}}, \code{\link{plot.core}}}
\examples{

\dontrun{
##-- Read kinesin alignment and structures
aln <- read.fasta(system.file("examples/kinesin_xray.fa",package="bio3d"))
pdb.path = system.file("examples/",package="bio3d")
pdbs <- read.fasta.pdb(aln, pdb.path = pdb.path, pdbext = ".ent")

## Find core
core <- core.find(pdbs,
                  #write.pdbs = TRUE,
                  verbose=TRUE)
}

##-- OR read previously saved kinesin data
data(kinesin)
attach(kinesin, warn.conflicts=FALSE)

col=rep("black", length(core$volume))
col[core$volume<2]="pink"; col[core$volume<1]="red"
plot(core, col=col)

inds.1 <- print(core, vol=1)
inds.0.5 <- print(core, vol=0.5)


}
\keyword{ utilities }

