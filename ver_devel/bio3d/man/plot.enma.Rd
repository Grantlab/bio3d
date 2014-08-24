\name{plot.enma}
\alias{plot.enma}
\title{ Plot eNMA Results }
\description{
  Produces a plot of atomic fluctuations of aligned normal modes.
}
\usage{
\method{plot}{enma}(x, 
        pdbs = NULL, conservation = NULL, variance = FALSE, 
        spread = FALSE, offset = 1,
        col = NULL, signif = FALSE,
        pcut = 0.005, qcut = 0.04,
        xlab = "Alignment Position",
        ylab=c("Fluctuations", "Fluct.variance", "Seq.conservation"),
        xlim = NULL, ylim = NULL,
        mar = c(4, 5, 2, 2), ...)
}
\arguments{
  \item{x}{ the results of ensemble NMA obtained with
    \code{\link{nma.pdbs}}. Alternatively, a matrix in the similar
    format as \code{enma$fluctuations} can be provided. }
  \item{pdbs}{ an object of class \sQuote{pdbs} in which the
    \sQuote{enma} object \code{x} was obtained from. If provided SSE
    data of the first structure of \code{pdbs} will drawn. }
  \item{conservation}{ logical, if TRUE sequence conservation is
    plotted. Alternatively, provide the conservation assement method
    (\sQuote{similarity}, \sQuote{identity}, \sQuote{entropy22}, or
    \sQuote{entropy10}). A numeric vector of residue conservation values
    are also allowed. }    
  \item{variance}{ logical, if TRUE fluctuation variance is plotted. }
  \item{spread}{ logical, if TRUE the fluctuation profiles are spread -
    i.e. not on top of each other. }
  \item{offset}{ numerical offset value in use when
    \sQuote{spread=TRUE}. }
  \item{col}{ a character vector of plotting colors. Used also to group
    fluctuation profiles when \sQuote{spread=TRUE}. NA values in col
    will omit the corresponding fluctuation profile in the plot. }
  \item{signif}{ logical, if TRUE significance of difference is plotted. }
  \item{pcut}{ P-value cutoff for the significance. }
  \item{qcut}{ Cutoff for the minimal difference of mean fluctuation to plot 
             the significance. }
  \item{xlab}{ a label for the x axis. }
  \item{ylab}{ labels for the y axes. }
  \item{mar}{ a numerical vector of the form c(bottom, left, top, right)
    which gives the number of lines of margin to be specified on the
    four sides of the plot.}
  \item{xlim}{ the x limits of the plot. }
  \item{ylim}{ the y limits of the plot. }
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
\dontrun{
ids <- c("1a70_A", "1czp_A", "1frd_A", "1fxi_A", "1iue_A", "1pfd_A")
raw.files <- get.pdb(ids, path = "raw_pdbs")
files <- pdbsplit(raw.files, ids, path = "raw_pdbs/split_chain")

## Sequence Alignement
pdbs <- pdbaln(files)

## Normal mode analysis on aligned data
all.modes <- nma.pdbs(pdbs, rm.gaps=TRUE)

## Plot fluctuations
plot.enma(all.modes, pdbs=pdbs, conservation=TRUE)

## group and spread fluctuation profiles
grps <- rep(NA, length(pdbs$id))
grps[c(2,3)]=1
grps[c(4,5)]=2

plot.enma(all.modes, pdbs=pdbs, col=grps, spread=TRUE)
}
}
\keyword{ hplot }
