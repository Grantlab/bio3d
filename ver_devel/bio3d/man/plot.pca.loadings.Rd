\name{plot.pca.loadings}
\alias{plot.pca.loadings}
\title{  Plot Residue Loadings along PC1 to PC3 }
\description{
   Plot residue loadings along PC1 to PC3 from a given xyz C-alpha matrix
   of \code{loadings}.
 }
\usage{
plot.pca.loadings(x, resnums = seq(1, (length(x[, 1])/3), 25), ...)
}
\arguments{
  \item{x}{ the results of principal component analysis obtained
    from \code{\link{pca.xyz}}, or just the loadings returned from
    \code{\link{pca.xyz}}. }
  \item{resnums}{ a numeric vector of residue numbers. }
  \item{\dots}{ extra plotting arguments. }  
}
\value{
  Called for its effect.
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Barry Grant }
\seealso{ \code{\link{pca.xyz}}, \code{\link{plot.pca}} }
\examples{
data(kinesin)
attach(kinesin)
plot.pca.loadings(pc.xray$U)
}
\keyword{ hplot }
