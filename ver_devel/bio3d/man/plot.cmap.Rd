\name{plot.cmap}
\alias{plot.cmap}
\title{ Plot Contact Matrix }
\description{
  Plot a contact matrix with optional secondary structure in the
  marginal regions.
}
\usage{
\method{plot}{cmap}(x, col=2, pch=16, main="Contact map", sub="",
           xlim=NULL, ylim=NULL, xlab = "Residue index", ylab = xlab, 
           axes=TRUE, ann=par("ann"), sse=NULL, sse.type="classic",
           sse.min.length=5, bot=TRUE, left=TRUE,
           helix.col="gray20", sheet.col="gray80", sse.border=FALSE,
           add=FALSE, ...)
}
\arguments{
  \item{x}{ a numeric matrix of residue contacts as obtained from
    function \code{cmap}. }
  \item{col}{ color code or name, see \code{par}. }
  \item{pch}{ plotting \sQuote{character}, i.e., symbol to use. This can
    either be a single character or an integer code for one of a set of
    graphics symbols. See \code{points}. }
  \item{main}{ a main title for the plot, see also \sQuote{title}. }
  \item{sub}{ a sub-title for the plot. }
  \item{xlim}{ the x limits (x1,x2) of the plot.  Note that x1 > x2 is
          allowed and leads to a reversed axis. }
  \item{ylim}{ the y limits of the plot. }
  \item{xlab}{ a label for the x axis, defaults to a description of \sQuote{x}. }
  \item{ylab}{ a label for the y axis, defaults to a description of \sQuote{y}. }
  \item{axes}{ a logical value indicating whether both axes should be drawn
    on the plot.  Use graphical parameter \sQuote{xaxt} or
    \sQuote{yaxt} to suppress just one of the axes. }
  \item{ann}{ a logical value indicating whether the default annotation
    (title and x and y axis labels) should appear on the plot. }
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
  \item{left}{ logical, if TRUE rectangles for each sse are drawn towards
    the left of the plotting region.  }
  \item{bot}{ logical, if TRUE rectangles for each sse are drawn towards
    the bottom of the plotting region. }
  \item{helix.col}{ The colors for rectangles representing alpha helices. }
  \item{sheet.col}{ The colors for rectangles representing beta strands. }
  \item{sse.border}{ The border color for all sse rectangles. }
  \item{add}{ logical, specifying if the contact map should be added to an
    already existing plot. Note that when \sQuote{TRUE} only points are
    plotted (no annotation). }
  \item{\dots}{ other graphical parameters. }
}
\details{
  This function is useful for plotting a residue-residue contact data
  for a given protein structure along with a schematic representation of
  major secondary structure elements. 

  Two forms of secondary structure annotation are available: so called
  \sQuote{classic} and \sQuote{fancy}. The former draws marginal
  rectangles and has been available within Bio3D from version 0.1. The
  later draws more \sQuote{fancy} (and distracting) 3D like helices and
  arrowed strands.
}
\value{
  Called for its effect.
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Lars Skjaerven, Barry Grant }
\note{ Be sure to check the correspondence of your \sQuote{sse} object with the
  \sQuote{x} values being plotted as no internal checks are performed.
}
\seealso{
  \code{\link{cmap}}, \code{\link{dm}},
  \code{\link{plot.dmat}}, 
  \code{\link{plot.default}}, \code{\link{plot.bio3d}},
  \code{\link{dssp}}, \code{\link{stride}}
}
\examples{

##- Read PDB file
pdb <- read.pdb( system.file("examples/1hel.pdb", package="bio3d") )

##- Calcualte contact map
cm <- cmap(pdb)

##- Plot contact map
plot.cmap(cm, sse=pdb)

##- Add to plot
plot.cmap(t(cm), col=3, pch=17, add=TRUE)

}
\keyword{ hplot }
