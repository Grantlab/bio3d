\name{plot.enma}
\alias{plot.enma}
\title{ Plot eNMA Results }
\description{
  Produces a plot of atomic fluctuations of aligned normal modes.
}
\usage{
\method{plot}{enma}(x, pdbs = NULL, entropy = FALSE, col = NULL,
        xlab = "Residue Position", ylab = "Fluctuations",
        mar = c(4, 5, 2, 2), ...)
}
\arguments{
  \item{x}{ the results of ensemble NMA obtained with
    \code{\link{nma.pdbs}}. }
  \item{pdbs}{ an object of class \sQuote{3dalign} in which the
    \sQuote{enma} object \code{x} was obtained from. If provided SSE
    data of the first structure of \code{pdbs} will drawn. }
  \item{entropy}{ logical, if TRUE entropy and fluctuation variance is
    plotted. }
  \item{col}{ a character vector of plotting colors.  }
  \item{xlab}{ a label for the x axis. }
  \item{ylab}{ a label for the y axis. }
  \item{mar}{ A numerical vector of the form c(bottom, left, top, right)
    which gives the number of lines of margin to be specified on the
    four sides of the plot.}
  \item{\dots}{ extra plotting arguments passed to \code{plot.bio3d}
    that effect the atomic fluctuations plot only. }
}
\details{
  \code{plot.enma} produces a fluctuation plot of aligned \code{nma}
  objects.
}
\value{
  Called for its effect.
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Lars Skjaerven, Barry Grant }
\seealso{ \code{\link{nma.pdbs}}, \code{\link{nma}},
  \code{\link{plot.bio3d}}, \code{\link{entropy}}.  }
\examples{
ids <- c("1a70_A", "1czp_A", "1frd_A", "1fxi_A", "1iue_A", "1pfd_A")
raw.files <- get.pdb(ids, path = "raw_pdbs")
files <- pdbsplit(raw.files, ids, path = "raw_pdbs/split_chain")

## Sequence Alignement
pdbs <- pdbaln(files)

## Normal mode analysis on aligned data
all.modes <- nma.pdbs(pdbs, rm.gaps=FALSE)

## Plot fluctuations
plot.enma(all.modes, pdbs)
}
\keyword{ hplot }
