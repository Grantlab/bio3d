\name{rmsd.filter}
\alias{rmsd.filter}
\title{ RMSD Filter }
\description{
  Identify and filter subsets of conformations at a given RMSD cutoff.
}
\usage{
rmsd.filter(xyz = NULL, rmsd.mat = NULL, cutoff = 0.5,
            fit = TRUE, verbose = TRUE)
}
\arguments{
  \item{xyz}{ a numeric matrix or list object containing multiple
    coordinates for pairwise comparison, such as that obtained from \code{\link{read.fasta.pdb}}. Not used if \code{rmsd.mat} is given. }
  \item{rmsd.mat}{  an optional matrix of RMSD values obtained from \code{\link{rmsd}}. }
  \item{cutoff}{ a numeric rmsd cutoff value. }
  \item{fit}{ logical, if TRUE coordinate superposition is performed
    prior to RMSD calculation. }
  \item{verbose}{ logical, if TRUE progress details are printed. }
}
\details{
  This function performs hierarchical cluster analysis of a given matrix of
  RMSD values \sQuote{rmsd.mat}, or an RMSD matrix calculated from a given
  coordinate matrix \sQuote{xyz}, to identify conformers that fall below a
  given RMSD cutoff value \sQuote{cutoff}.
}
\value{
  Returns a list object with components:
  \item{ind}{indices of the conformers (rows) below the cutoff value.}
  \item{tree}{an object of class \code{"hclust"}, which describes the
    tree produced by the clustering process. }
  \item{rmsd.mat}{a numeric matrix with all pairwise RMSD values.}
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Barry Grant }
\seealso{ \code{\link{rmsd}}, \code{\link{read.pdb}},
  \code{\link{read.fasta.pdb}}, \code{\link{read.dcd}}
}
\examples{
\dontrun{
data(kinesin)
attach(kinesin, warn.conflicts=FALSE)
k <- rmsd.filter(xyz=pdbs,cutoff=0.5)
aln$id[k$ind]
plot(k$tree, ylab="RMSD")
abline(h=0.5, col="gray")
}
}
\keyword{ utilities }
