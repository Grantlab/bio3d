\name{plot.enma}
\alias{plot.enma}
\title{ Plot eNMA Results }
\description{
  Produces a plot of atomic fluctuations of aligned normal modes.
}
\usage{
\method{plot}{enma}(x,
        pdbs = NULL,
        xlab = NULL,
        ylab="Fluctuations", ...)
}
\arguments{
  \item{x}{ the results of ensemble NMA obtained with
    \code{\link{nma.pdbs}}. Alternatively, a matrix in the similar
    format as \code{enma$fluctuations} can be provided. }
  \item{pdbs}{ an object of class \sQuote{pdbs} in which the
    \sQuote{enma} object \code{x} was obtained from. If provided SSE
    data of the first structure of \code{pdbs} will drawn. }
  \item{xlab}{ a label for the x axis. }
  \item{ylab}{ labels for the y axes. }
  \item{\dots}{ extra plotting arguments passed to \code{plot.fluct}
    that effect the atomic fluctuations plot only. }
}
\details{
  \code{plot.enma} produces a fluctuation plot of aligned \code{nma}
  objects. If corresponding \code{pdbs} object is provided the plot
  contains SSE annotation and appropriate resiude index numbering.
}
\value{
  Called for its effect.
}
\references{
  Skjaerven, L. et al. (2014) \emph{BMC Bioinformatics} \bold{15}, 399.
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Lars Skjaerven, Barry Grant }
\seealso{ \code{\link{nma}},
  \code{\link{plotb3}}, \code{\link{plot.fluct}} }
\examples{
\dontrun{
ids <- c("1a70_A", "1czp_A", "1frd_A", "1fxi_A", "1iue_A", "1pfd_A")
raw.files <- get.pdb(ids, path = "raw_pdbs")
files <- pdbsplit(raw.files, ids, path = "raw_pdbs/split_chain")

## Sequence/structure alignement
pdbs <- pdbaln(files)

## Normal mode analysis on aligned data
modes <- nma(pdbs)

## Plot fluctuations
plot(modes, pdbs=pdbs)

## Group and spread fluctuation profiles
hc <- hclust(as.dist(1-modes$rmsip))
col <- cutree(hc, k=2)
plot(modes, pdbs=pdbs, col=col, spread=TRUE)
}
}
\keyword{ hplot }
