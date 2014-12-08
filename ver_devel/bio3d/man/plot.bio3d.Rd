\name{plot.bio3d}
\alias{plot.bio3d}
\title{ Plots with marginal SSE annotation }
\description{
  Draw a standard scatter plot with optional secondary structure in the
  marginal regions.
}
\usage{
  plot.bio3d(x, resno = NULL, rm.gaps = FALSE, type = "h",
             main = "", sub = "",
             xlim = NULL, ylim = NULL, ylim2zero = TRUE,
             xlab = "Residue", ylab = NULL, 
             axes = TRUE, ann = par("ann"), col = par("col"),
             sse = NULL, sse.type="classic", sse.min.length=5,
             top = TRUE, bot = TRUE, 
             helix.col = "gray20", sheet.col = "gray80",
             sse.border = FALSE, ...)
}
\arguments{
  \item{x}{ a numeric vector of values to be plotted.  Any reasonable way of
    defining these plot coordinates is acceptable.  See the function
    \sQuote{xy.coords} for details. }
  \item{resno}{ an optional vector with length equal to that of
    \sQuote{x} that will be used to annotate the xaxis. This is typically
    a vector of residue numbers. If NULL residue positions from 1 to the
    length of \sQuote{x} will be used. See examples below. }
  \item{rm.gaps}{ logical, if TRUE gaps in \code{x}, indicated by NA values,
    will be removed from plot. }
  \item{type}{ one-character string giving the type of plot desired.  The
    following values are possible, (for details, see \sQuote{plot}): \sQuote{p}
    for points, \sQuote{l} for lines, \sQuote{o} for over-plotted points and
    lines, \sQuote{b}, \sQuote{c}) for points joined by lines,
    \sQuote{s} and \sQuote{S} for stair steps and \sQuote{h} for
    histogram-like vertical lines.  Finally, \sQuote{n} does not produce
    any points or lines. }
  \item{main}{ a main title for the plot, see also \sQuote{title}. }
  \item{sub}{ a sub-title for the plot. }
  \item{xlim}{ the x limits (x1,x2) of the plot.  Note that x1 > x2 is
          allowed and leads to a reversed axis. }
  \item{ylim}{ the y limits of the plot. }
  \item{ylim2zero}{ logical, if TRUE the y-limits are forced to start at zero. }
  \item{xlab}{ a label for the x axis, defaults to a description of \sQuote{x}. }
  \item{ylab}{ a label for the y axis, defaults to a description of \sQuote{y}. }
  \item{axes}{ a logical value indicating whether both axes should be drawn
    on the plot.  Use graphical parameter \sQuote{xaxt} or
    \sQuote{yaxt} to suppress just one of the axes. }
  \item{ann}{ a logical value indicating whether the default annotation
    (title and x and y axis labels) should appear on the plot. }
  \item{col}{ The colors for lines and points.  Multiple colors can be
    specified so that each point is given its own color.  If there are
    fewer colors than points they are recycled in the standard fashion.
    Lines are plotted in the first color specified. }
  \item{sse}{ secondary structure object as returned from
    \code{\link{dssp}}, \code{\link{stride}} or in certain cases
    \code{\link{read.pdb}}. }
  \item{sse.type}{ single element character vector that determines the
    type of secondary structure annotation drawn. The following values
    are possible, \sQuote{classic} and \sQuote{fancy}. See details and 
    examples below. }
  \item{sse.min.length}{ a single numeric value giving the length below
    which secondary structure elements will not be drawn. This is useful
    for the exclusion of short helix and strand regions that can often crowd
    these forms of plots. }
  \item{top}{ logical, if TRUE rectangles for each sse are drawn towards
    the top of the plotting region.  }
  \item{bot}{ logical, if TRUE rectangles for each sse are drawn towards
    the bottom of the plotting region. }
  \item{helix.col}{ The colors for rectangles representing alpha helices. }
  \item{sheet.col}{ The colors for rectangles representing beta strands. }
  \item{sse.border}{ The border color for all sse rectangles. }
  \item{\dots}{ other graphical parameters. }
}
\details{
  This function is useful for plotting per-residue numeric vectors for a
  given protein structure (e.g. results from RMSF, PCA, NMA etc.) along
  with a schematic representation of major secondary structure elements.

  Two forms of secondary structure annotation are available: so called
  \sQuote{classic} and \sQuote{fancy}. The former draws marginal
  rectangles and has been available within Bio3D from version 0.1. The
  later draws more \sQuote{fancy} (and distracting) 3D like helices and
  arrowed strands.
  
  See the functions \sQuote{plot.default}, \code{\link{dssp}} and \code{\link{stride}} 
  for further details.
}
\value{
  Called for its effect.
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Barry Grant }
\note{ Be sure to check the correspondence of your \sQuote{sse} object with the
  \sQuote{x} values being plotted as no internal checks are performed.
}
\seealso{ \code{\link{plot.default}}, \code{\link{dssp}}, \code{\link{stride}} }
\examples{
\donttest{
## Plot of B-factor values along with secondary structure from PDB
pdb <- read.pdb( "1bg2" )
bfac <- pdb$atom[pdb$calpha,"b"]
plot.bio3d(bfac, sse=pdb, ylab="B-factor", col="gray")
points(bfac, typ="l")
}

\dontrun{
## Use PDB residue numbers and include short secondary structure elements
plot.bio3d(pdb$atom[pdb$calpha,"b"], sse=pdb, resno=pdb, ylab="B-factor",
  typ="l", lwd=1.5, col="blue", sse.min.length=0)


## Calculate secondary structure using stride() or dssp()
#sse <- stride(pdb)
sse <- dssp(pdb)

## Plot of B-factor values along with calculated secondary structure
plot.bio3d(pdb$atom[pdb$calpha,"b"], sse=sse, ylab="B-factor", typ="l",
col="blue", lwd=2)
}

## Plot 'aligned' data respecting gap positions
attach(transducin)

pdb = read.pdb("1tnd") ## Reference PDB see: pdbs$id[1]
pdb = trim.pdb(pdb, inds=atom.select(pdb, chain="A"))

## Plot of B-factor values with gaps
plot.bio3d(pdbs$b, resno=pdb, sse=pdb, ylab="B-factor")

## Plot of B-factor values after removing all gaps 
plot.bio3d(pdbs$b, rm.gaps=TRUE, resno = pdb, sse=pdb, ylab="B-factor")

detach(transducin)

## Fancy secondary structure elements
##plot.bio3d(pdb$atom[pdb$calpha,"b"], sse=pdb, ssetype="fancy")
## Currently not implemented
}
\keyword{ hplot }
