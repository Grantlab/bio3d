\name{rmsf}
\alias{rmsf}
\title{ Atomic RMS Fluctuations }
\description{
  Calculate atomic root mean squared fluctuations.
}
\usage{
rmsf(xyz, average=FALSE)
}
\arguments{
  \item{xyz}{ numeric matrix of coordinates with each row corresponding
    to an individual conformer. }
  \item{average}{ logical, if TRUE averaged over atoms. }
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
attach(transducin)

# Ignore Gaps
gaps <- gap.inspect(pdbs$ali)

r <- rmsf(pdbs$xyz)
plot(r[gaps$f.inds], typ="h", ylab="RMSF (A)")

detach(transducin)
}
\keyword{ utilities }
