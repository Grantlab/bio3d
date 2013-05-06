\name{plot.pca}
\alias{plot.pca}
\alias{plot.pca.score}
\alias{plot.pca.scree}
\title{ Plot PCA Results }
\description{
  Produces a z-score plot (conformer plot) and an eigen spectrum plot
  (scree plot).
}
\usage{
plot.pca(x, pch = 16, col = par("col"), cex=0.8, mar=c(4, 4, 1, 1),...)
plot.pca.scree(x, y = NULL, type = "o", pch = 18,
         main = "", sub = "", xlim = c(0, 20), ylim = NULL,
         ylab = "Proporton of Variance (\%)",
         xlab = "Eigenvalue Rank", axes = TRUE, ann = par("ann"),
         col = par("col"), lab = TRUE, ...)
plot.pca.score(x, start.row = 1, end.row = nrow(x), col=rainbow(nrow(x)), lab = "", ...)

}
\arguments{
  \item{x}{ the results of principal component analysis obtained with
    \code{\link{pca.xyz}}. }
  \item{pch}{ a vector of plotting characters or symbols: see \sQuote{points}. }
  \item{col}{ a character vector of plotting colors.  }
  \item{cex}{ a numerical single element vector giving the amount by
    which plotting text and symbols should be magnified relative to the
    default. }
  \item{mar}{ A numerical vector of the form c(bottom, left, top, right)
    which gives the number of lines of margin to be specified on the
    four sides of the plot.}
  \item{start.row}{ row index of the first conformer to colour. }
  \item{end.row}{ row index of the last conformer to colour. }
  \item{lab}{ a character vector of plot labels. }
  \item{y}{ the y coordinates for the scree plot. }
  \item{type}{ one-character string giving the type of plot desired.}
  \item{main}{ a main title for the plot, see also 'title'.}
  \item{sub}{ a sub-title for the plot.}
  \item{xlim}{ the x limits of the plot. }
  \item{ylim}{ the y limits of the plot.}
  \item{ylab}{ a label for the y axis.}
  \item{xlab}{ a label for the x axis.}
  \item{axes}{ a logical value indicating whether both axes should be drawn.}
  \item{ann}{ a logical value indicating whether the default annotation
    (title and x and y axis labels) should appear on the plot. }
  \item{\dots}{ extra plotting arguments. }
}
\details{
  \code{plot.pca} is a wrapper calling both \code{plot.pca.score} and
  \code{plot.pca.scree} resulting in a 2x2 plot with three score plots
  and one scree plot.
}
\value{
  Called for its effect.
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Barry Grant }
\seealso{ \code{\link{pca.xyz}}, \code{\link{plot.bio3d}} }
\examples{
data(kinesin)
attach(kinesin, warn.conflicts=FALSE)

plot(pc.xray)

## change colors
plot(pc.xray,col=rainbow(length(pc.xray$z[,1])), cex=1.5)

## add labels

## color by seq identity


}
\keyword{ hplot }
