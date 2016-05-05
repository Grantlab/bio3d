\name{plot.fluct}
\alias{plot.fluct}
\title{ Plot Fluctuations }
\description{
  Produces a plot of atomic fluctuations obtained from ensemble normal mode analysis or
  molecular dynamics simulations.
}
\usage{
\method{plot}{fluct}(x,
        col = NULL, label = rownames(x), signif = FALSE,
        p.cutoff = 0.005, q.cutoff = 0.04,
        s.cutoff = 5, n.cutoff = 2, mean = FALSE, polygon = FALSE,
        spread = FALSE, offset = 1,
        ncore = NULL, ...)
}
\arguments{
  \item{x}{ a numeric vector or matrix containing atomic fluctuation data obtained
    from e.g. \code{\link{nma.pdbs}} or \code{\link{rmsf}}. }
  \item{col}{ a character vector of plotting colors. Used also to group
    fluctuation profiles. NA values in col will omit the corresponding fluctuation
    profile in the plot. }
  \item{label}{ a character vector of plotting labels with length matching
    \code{nrow(x)}. If \code{mean=TRUE}, the length of \code{label} can be equal to
    the number of categories indicated by \code{col}. }
  \item{signif}{ logical, if TRUE significance of fluctuation difference is calculated
             and annotated for each atomic position. }
  \item{p.cutoff}{ Cutoff of p-value to define significance. }
  \item{q.cutoff}{ Cutoff of the mean fluctuation difference to define significance. }
  \item{s.cutoff}{ Cutoff of sample size in each group to calculate the significance. }
  \item{n.cutoff}{ Cutoff of consecutive residue positions with significant fluctuation
             difference. If the actual number is less than the cutoff,
             correponding postions will not be annotated. }
  \item{mean}{ logical, if TRUE plot mean fluctuations of each group. Significance is
             still calculated with the original data. }
  \item{polygon}{ logical, if TRUE a nicer plot with area under the line for the first
             row of \code{x} are filled with polygons. }
  \item{ncore }{ number of CPU cores used to do the calculation. By default
             (\code{ncore=NULL}), use all available CPU cores. The argument is only
             used when \code{signif=TRUE}. }
  \item{spread}{ logical, if TRUE the fluctuation profiles are spread -
    i.e. not on top of each other. }
  \item{offset}{ numerical offset value in use when
    \sQuote{spread=TRUE}. }
  \item{\dots}{ extra plotting arguments passed to \code{plot.bio3d}. }
}
\details{
  The significance calculation is performed when \code{signif=TRUE} and there are at least
  two groups with sample size larger than or equal to \code{s.cutoff}. A "two-sided"
  student's t-test is performed for each atomic position (each
  column of \code{x}).  If \code{x} contains gaps, indicated by \code{NA}s,
  only non-gapped positions are considered. The position is considered significant if both
  \code{p-value <= p.cutoff} and the mean value difference of the two groups, q, satisfies
  \code{q >= q.cutoff}. If more than two groups are available, every pair of groups are
  subjected to the t-test calculation and the minimal p-value along with the q-value
  for the corresponding pair are used for the significance evaluation.
}
\value{
  If significance is calculated, return a vector indicating significant positions.
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Xin-Qiu Yao, Lars Skjaerven, Barry Grant }
\seealso{ \code{\link{plot.bio3d}}, \code{\link{rmsf}}, \code{\link{nma.pdbs}},
          \code{\link[stats:t.test]{t.test}}, \code{\link[graphics:polygon]{polygon}}. }
\examples{
\dontrun{
## load transducin example data
attach(transducin)

## subset of pdbs to analyze
inds = c(1:5, 16:20)
pdbs <- trim(pdbs, row.inds=inds)
gaps.res = gap.inspect(pdbs$ali)

## reference RESNO and SSE for axis annotations
resno <- pdbs$resno[1, gaps.res$f.inds]
sse <- pdbs$sse[1, gaps.res$f.inds]

## eNMA calculation and obtain modes of motion including atomic fluctuations
modes <- nma(pdbs, ncore=NULL)
x = modes$fluctuation

## simple line plot with SSE annotation
plot.fluct(x, sse=sse, resno=resno)

## group data by specifying colors of each fluctuation line; same color indicates
## same group. Also do significance calculation and annotation
col = c(rep('red', 5), rep('blue', 5))
plot.fluct(x, col=col, signif=TRUE, sse=sse, resno=resno)

## spread lines
plot.fluct(x, col=col, signif=TRUE, sse=sse, resno=resno, typ='l', spread=TRUE)

## show only line of mean values for each group. 
## Nicer plot with area shaded for the first group.
plot.fluct(x, col=col, signif=TRUE, sse=sse, resno=resno, mean=TRUE, 
           polygon=TRUE, label=c('GTP', 'GDI'))

detach(transducin)
}
}
\keyword{ hplot }
