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
  RMSF is an often used measure of conformational variance. It is calculated by
  \deqn{f_i=\sqrt{\frac{1}{M-1}\sum_j \|r_i^j-r_i^0\|^2}}{%
        f_i=sqrt((1/(M-1))Sum_j(||r_i^j-r_i^0||^2))},
  where \eqn{f_i}{f_i} is the RMSF value for the ith atom, M the total number of frames 
  (total number of rows of \code{xyz}), \eqn{r_i^j}{r_i^j} the positional vector of the
  ith atom in the jth frame, and \eqn{r_i^0}{r_i^0} the mean position of ith atom. 
  ||r|| denotes the Euclidean norm of the vector r.
}
\value{
  Returns a numeric vector of RMSF values. If \code{average=TRUE} a single numeric value
  representing the averaged RMSF value over all atoms will be returned. 
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
