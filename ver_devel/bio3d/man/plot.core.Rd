\name{plot.core}
\alias{plot.core}
\title{ Plot Core Fitting Progress }
\description{
  Plots the total ellipsoid volume of core positions versus core size at
  each iteration of the core finding process.
}
\usage{
\method{plot}{core}(x, y = NULL, type = "h", main = "", sub = "", xlim = NULL, ylim = NULL, xlab = "Core Size (Number of Residues)", ylab = "Total Ellipsoid Volume (Angstrom^3)", axes = TRUE, ann = par("ann"), col = par("col"), ...)
}
\arguments{
  \item{x}{  a list object obtained with the function
    \code{\link{core.find}} from which the \sQuote{volume}
    component is taken as the x coordinates for the plot. }
  \item{y}{ the y coordinates for the plot. }
  \item{type}{ one-character string giving the type of plot desired. }
  \item{main}{ a main title for the plot, see also \sQuote{title}. }
  \item{sub}{ a sub-title for the plot. }
  \item{xlim}{ the x limits of the plot. }
  \item{ylim}{ the y limits of the plot. }
  \item{xlab}{ a label for the x axis. }
  \item{ylab}{ a label for the y axis. }
  \item{axes}{ a logical value indicating whether both axes should be drawn. }
  \item{ann}{ a logical value indicating whether the default annotation
    (title and x and y axis labels) should appear on the plot. }
  \item{col}{ The colors for lines and points.  Multiple colours can be
    specified so that each point is given its own color.  If there are
    fewer colors than points they are recycled in the standard fashion. }
  \item{\dots}{ extra plotting arguments. }
}
\value{
  Called for its effect.
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Barry Grant }
\note{
  The produced plot can be useful for deciding on the core/non-core boundary.
}
\seealso{ \code{\link{core.find}}, \code{\link{print.core}}}
\examples{
\dontrun{

##-- Generate a small kinesin alignment and read corresponding structures
pdbfiles <- get.pdb(c("1bg2","2ncd","1i6i","1i5s"), URLonly=TRUE)
pdbs <- pdbaln(pdbfiles)

##-- Find 'core' positions
core <- core.find(pdbs)
plot(core)

##-- Fit on these relatively invarient subset of positions 
core.inds <- print(core)
xyz <- pdbfit(pdbs, core.inds, outpath="corefit_structures/")

##-- Compare to fitting on all equivalent positions
xyz2 <- pdbfit(pdbs)

## Note that overall RMSD will be higher but RMSF will
##  be lower in core regions, which may equate to a
##  'better fit' for certain applications
gaps <- gap.inspect(pdbs$xyz)
rmsd(xyz[,gaps$f.inds])
rmsd(xyz2[,gaps$f.inds])

plot(rmsf(xyz[,gaps$f.inds]), typ="l", col="blue", ylim=c(0,9))
points(rmsf(xyz2[,gaps$f.inds]), typ="l", col="red")

}

}
\keyword{ hplot }
