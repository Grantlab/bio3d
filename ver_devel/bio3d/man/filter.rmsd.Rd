\name{filter.rmsd}
\alias{filter.rmsd}
\title{ RMSD Filter }
\description{
  Identify and filter subsets of conformations at a given RMSD cutoff.
}
\usage{
filter.rmsd(xyz = NULL, rmsd.mat = NULL, cutoff = 0.5,
            fit = TRUE, verbose = TRUE, inds = NULL, method = "complete",
            \dots)
}
\arguments{
  \item{xyz}{ a numeric matrix or list object containing multiple
    coordinates for pairwise comparison, such as that obtained from
    \code{\link{read.fasta.pdb}}. Not used if \code{rmsd.mat} is
    given. }
  \item{rmsd.mat}{  an optional matrix of RMSD values obtained from
    \code{\link{rmsd}}. } 
  \item{cutoff}{ a numeric rmsd cutoff value. }
  \item{fit}{ logical, if TRUE coordinate superposition is performed
    prior to RMSD calculation. }
  \item{verbose}{ logical, if TRUE progress details are printed. }
  \item{inds}{ a vector of indices that selects the elements of
    \code{xyz} upon which the calculation should be based. By default, 
    all the non-gap sites in \code{xyz}. }
  \item{method}{ the agglomeration method to be used. See function 
    \code{\link{hclust}} for more information. }
  \item{\dots}{ additional arguments passed to and from functions. } 
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
\seealso{
  \code{\link{rmsd}}, \code{\link{read.pdb}},
  \code{\link{read.fasta.pdb}}, \code{\link{read.dcd}}
}
\examples{
\dontrun{
data(kinesin)
attach(kinesin, warn.conflicts=FALSE)
k <- filter.rmsd(xyz=pdbs,cutoff=0.5)
pdbs$id[k$ind]
plot(k$tree, ylab="RMSD")
abline(h=0.5, col="gray")

detach(kinesin)
}
}
\keyword{ utilities }
