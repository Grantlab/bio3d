\name{rmsf}
\alias{rmsf}
\title{ Atomic RMS Fluctuations }
\description{
  Calculate atomic root mean squared fluctuations.
}
\usage{
rmsf(xyz, grpby=NULL, average=FALSE)
}
\arguments{
  \item{xyz}{ numeric matrix of coordinates with each row corresponding
    to an individual conformer. }
  \item{grpby}{ a vector counting connective duplicated elements that
    indicate the elements of 'xyz' that should be considered as a
    group (e.g. atoms from a particular residue). If provided a
    'pdb' object, grouping is automatically set by amino acid residues. }
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

\dontrun{

pdb <- read.pdb("1d1d", multi=TRUE)
xyz <- pdb$xyz

# superimpose trajectory
xyz <- fit.xyz(xyz[1, ], xyz)

# select mainchain atoms
sele <- atom.select(pdb, elety=c("CA", "C", "N", "O"))

# residue numbers to group by
resno <- pdb$atom$resno[sele$atom]

# mean rmsf value of mainchain atoms of each residue
r <- rmsf(xyz[, sele$xyz], grpby=resno)
plot.bio3d(r, resno=pdb, sse=pdb, ylab="RMSF (A)")

}

}
\keyword{ utilities }
