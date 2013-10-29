\name{plot.enma}
\alias{plot.enma}
\title{ Plot eNMA Results }
\description{
  Produces a plot of atomic fluctuations of aligned normal modes.
}
\usage{
\method{plot}{enma}(x, xlab = "Residue Position", ylab = "Fluctuations",
                col = NULL, cex = 0.8, mar = c(6, 5, 2, 2), ...)
}
\arguments{
  \item{x}{ the results of ensemble NMA obtained with
    \code{\link{nma.pdbs}}. }
  \item{xlab}{ a label for the x axis, defaults to a description of \sQuote{x}. }
  \item{ylab}{ a label for the y axis, defaults to a description of \sQuote{y}. }
  \item{col}{ a character vector of plotting colors.  }
  \item{cex}{ a numerical single element vector giving the amount by
    which plotting text and symbols should be magnified relative to the
    default. }
  \item{mar}{ A numerical vector of the form c(bottom, left, top, right)
    which gives the number of lines of margin to be specified on the
    four sides of the plot.}
  \item{\dots}{ extra plotting arguments passed to \code{plot.bio3d}
    that effect the atomic fluctuations plot only. }
}
\details{
  \code{plot.enma} produces a fluctuation plot of aligned normal mode
  vectors. 
}
\value{
  Called for its effect.
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Lars Skjaerven, Barry Grant }
\seealso{ \code{\link{nma.pdbs}}, \code{\link{nma}},
  \code{\link{plot.bio3d}}.  }
\examples{
ids <- c("1a70_A", "1czp_A", "1frd_A", "1fxi_A", "1iue_A", "1pfd_A")
raw.files <- get.pdb(ids, path = "raw_pdbs")
files <- pdbsplit(raw.files, ids, path = "raw_pdbs/split_chain")

## Sequence Alignement
pdbs <- pdbaln(files)

## Normal mode analysis on aligned data
all.modes <- nma.pdbs(pdbs, rm.gaps=FALSE)

## Plot fluctuations
plot.enma(all.modes)

}
\keyword{ hplot }
