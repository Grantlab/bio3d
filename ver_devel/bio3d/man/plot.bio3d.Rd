\name{plot.bio3d}
\alias{plot.bio3d}
\title{ Plots with marginal SSE annotation }
\description{
  Draw a standard scatter plot with optional secondary structure in the
  marginal regions.
}
\usage{
\method{plot}{bio3d}(x, y = NULL, type = "h", main = "", sub = "", xlim = NULL, ylim = NULL, 
                 ylim2zero = TRUE, xlab = NULL, ylab = NULL, axes = TRUE, 
                 ann = par("ann"), col = par("col"), sse = NULL, top = TRUE, bot = TRUE, 
                 helix.col = "gray20", sheet.col = "gray80", sse.border = FALSE, ...)
}
\arguments{
  \item{x}{ the x coordinates for the plot.  Any reasonable way of
    defining the coordinates is acceptable.  See the function
    \sQuote{xy.coords} for details. }
  \item{y}{ the y coordinates for the plot, see above. }
  \item{type}{ one-character string giving the type of plot desired.  The
    following values are possible, (for details, see \sQuote{plot}): \sQuote{p}
    for points, \sQuote{l} for lines, \sQuote{o} for overplotted points and
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
  \item{col}{ The colors for lines and points.  Multiple colours can be
    specified so that each point is given its own color.  If there are
    fewer colors than points they are recycled in the standard fashion.
    Lines are plotted in the first colour specified. }
  \item{sse}{ secondary structure object as returned from
    \code{\link{dssp}} or \code{\link{stride}}. }
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
\dontrun{
##-- SSE setup --#
pdb <- read.pdb( "1bg2" )
sse <- stride(pdb, resno=FALSE)
}

##-- Or load previous dataset
data(kinesin)
attach(kinesin, warn.conflicts=FALSE)

## Plot of B-factor values along with secondary structure
plot.bio3d(pdb$atom[pdb$calpha,"b"], sse=sse, ylab="B-factor")

plot.bio3d(pdb$atom[pdb$calpha,"b"], sse=sse, ylab="B-factor", typ="l",
col="blue", lwd=2)

}
\keyword{ hplot }
