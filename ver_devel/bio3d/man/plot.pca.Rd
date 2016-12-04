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
\method{plot}{pca}(x, pc.axes=NULL, pch=16, col=par("col"), cex=0.8, mar=c(4, 4, 1, 1),...)

\method{plot}{pca.scree}(x, y = NULL, type = "o", pch = 18,
         main = "", sub = "", xlim = c(0, 20), ylim = NULL,
         ylab = "Proportion of Variance (\%)",
         xlab = "Eigenvalue Rank", axes = TRUE, ann = par("ann"),
         col = par("col"), lab = TRUE, ...)
\method{plot}{pca.score}(x, inds=NULL, col=rainbow(nrow(x)), lab = "", ...)

}
\arguments{
  \item{x}{ the results of principal component analysis obtained with
    \code{\link{pca.xyz}}. }
  \item{pc.axes}{ an optional numeric vector of length two specifying the principal components to be plotted.  A NULL value will result in an overview plot of the first three PCs and a scree plot. See examples. }    
  \item{pch}{ a vector of plotting characters or symbols: see \sQuote{points}. }
  \item{col}{ a character vector of plotting colors.  }
  \item{cex}{ a numerical single element vector giving the amount by
    which plotting text and symbols should be magnified relative to the
    default. }
  \item{mar}{ A numerical vector of the form c(bottom, left, top, right)
    which gives the number of lines of margin to be specified on the
    four sides of the plot.}
  \item{inds}{ row indices of the conformers to label. }
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
  Produces a plot of PCA results in the active graphics device and invisibly 
  returns the plotted \sQuote{z} coordinates along the requested \sQuote{pc.axes}. 
  See examples section where these coordinates are used to identify plotted points.
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Barry Grant }
\seealso{ \code{\link{pca.xyz}}, \code{\link{plot.bio3d}} }
\examples{
attach(transducin)

pc.xray <- pca(pdbs$xyz, rm.gaps=TRUE)
plot(pc.xray)

## Color plot by nucleotide state
vcolors <- annotation[, "color"]
plot(pc.xray, col=vcolors)

## Focus on a single plot of PC1 vs PC2
x <- plot(pc.xray, pc.axes=1:2, col=vcolors)

## Identify points interactively with mouse clicks
#identify(x, labels=basename.pdb(pdbs$id))

## Add labels to select points
inds <- c(1,10,37)
text(x[inds,], labels=basename.pdb(pdbs$id[inds]), col="blue")

## Alternative labeling method
#labs <- rownames(annotation)
#inds <- c(2,7)
#plot.pca.score(pc.xray, inds=inds, col=vcolors, lab=labs)

## color by seq identity groupings
#ide <- seqidentity(pdbs$ali)
#hc <- hclust(as.dist(1-ide))
#grps <- cutree(hc, h=0.2)
#vcolors <- rainbow(max(grps))[grps]
#plot.pca.score(pc.xray, inds=inds, col=vcolors, lab=labs)

detach(transducin)
}
\keyword{ hplot }
