\name{print.core}
\alias{print.core}
\title{ Printing Core Positions and Returning Indices }
\description{
  Print method for core.find objects.
}
\usage{
print.core(x, vol = NULL, ...)
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
##-- Generate a small kinesin alignment and read corresponding structures
pdbfiles <- get.pdb(c("1bg2","2ncd","1i6i","1i5s"), URLonly=TRUE)
pdbs <- pdbaln(pdbfiles)

##-- Find 'core' positions
core <- core.find(pdbs)
plot(core)

##-- Fit on these relatively invarient subset of positions 
core.inds <- print(core, vol=0.5)

print(core, vol=0.7)
print(core, vol=1.0)


}
\keyword{ utilities }

