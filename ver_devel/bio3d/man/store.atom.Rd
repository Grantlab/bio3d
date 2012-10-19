\name{store.atom}
\alias{store.atom}
\title{ Store all-atom data from a PDB object }
\description{
  Not intended for public usage
}
\usage{
store.atom(pdb)
}
\arguments{
  \item{pdb}{ A pdb object as obtained from read.pdb }
}
\details{
  This function was requested by a user and has not been extensively
  tested.  Hence it is not yet recommended for public usage.
}
\value{
  Returns a matrix of all-atom data
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Barry Grant }
\note{ This function is still in development and is NOT part of the
  offical bio3d package
}
\seealso{ \code{\link{read.fasta.pdb}} }
\examples{
\dontrun{
sel <- c("N", "CA", "C", "O", "CB", "*G", "*D",  "*E", "*Z")
pdbs <- read.all(aln, pdb.path = pdb.path, pdbext = ".ent", sel=sel)
core <- core.find(pdbs, rm.island=TRUE)

##-- Core indices
atm <- rep( rep(sel,each=3), ncol(aln$ali))
ca.ind  <- which(atm == "CA")
core.ind <- c( matrix(ca.ind, nrow=3)[,core$c0.5A.atom] )

##-- Core Fit
xyz <- fit.xyz(fixed=tmp[1,],
             mobile = tmp,
             fixed.inds  = core.ind,
             mobile.inds = core.ind)


##-- Examine all-atom data
colnames(tmp) = atm

gap.ind <- unique(which(is.na(tmp),arr.ind=T)[,2])
nogap.ind    <- c(1:ncol(tmp))[-gap.ind]
#ca.nogap.ind <- intersect(ca.ind, nogap.ind)


##-- PCA 
xray.pca <- pca.xyz( xyz[,nogap.ind] )

plot(xray.pca)
}
}
\keyword{ utilities }
