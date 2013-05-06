\name{rmsf}
\alias{rmsf}
\title{ Atomic RMS Fluctuations }
\description{
  Calculate atomic root mean squared fluctuations.
}
\usage{
rmsf(xyz)
}
\arguments{
  \item{xyz}{ numeric matrix of coordinates with each row corresponding
    to an individual conformer. }
}
\details{
  An often used measure of conformational variance.
}
\value{
  Returns a numeric vector of RMSF values.
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Barry Grant }
\seealso{ \code{\link{read.dcd}}, \code{\link{fit.xyz}},
  \code{\link{read.fasta.pdb}} }
\examples{
data(kinesin)
attach(kinesin, warn.conflicts=FALSE)

r <- rmsf(xyz)
plot( r[aln$ali["d1bg2__",]!="-"], typ="h", ylab="RMSF (A)")
}
\keyword{ utilities }
