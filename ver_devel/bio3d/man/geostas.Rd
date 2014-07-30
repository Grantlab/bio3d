\name{geostas}
\alias{geostas}
\alias{amsm.xyz}
\title{ GeoStaS Domain Finder }
\description{
  Identifies geometrically stable domains in biomolecules
}
\usage{
geostas(xyz, amsm = NULL, k = 3, method = "pairwise", fit = TRUE, \dots)

amsm.xyz(xyz, ncore = NULL)
}
\arguments{
  \item{xyz}{ numeric matrix of xyz coordinates as obtained e.g. by
    \code{\link{read.ncdf}}, \code{\link{read.dcd}}, or
    \code{\link{mktrj}}. }
  \item{amsm}{ a numeric matrix as obtained by
    \code{\link{amsm.xyz}} (convenient e.g. for re-doing only the
    clustering analysis of the \sQuote{AMSM} matrix). }
  \item{k }{ an integer scalar or vector with the desired number of
    groups. }
  \item{method }{ character string specifing the approach for clustering
    the atomic movement similarity matrix (pariwise or columnwise). }
  \item{fit}{ logical, if TRUE coordinate superposition on identified
    core atoms is performed prior to the calculation of the AMS
    matrix. } 
  \item{ncore }{ number of CPU cores used to do the calculation.
    \code{ncore>1} requires package \sQuote{parallel} installed. }
  \item{...}{ additional arguments to \code{\link{amsm.xyz}}. }
}
\details{
  This function attempts to identify rigid domains in a protein (or
  nucleic acid) structure based on an structural ensemble, e.g. obtained
  from NMR experiments, molecular dynamics simulations, or normal mode
  analysis.
  
  The algorithm is based on a geometric approach for comparing pairwise
  traces of atomic motion and the search for their best superposition
  using a quaternion representation of rotation. The result is stored in
  a NxN atomic movement similarity matrix (AMSM) describing the
  correspondence between all pairs of atom motion. Rigid domains are
  obtained by clustering the elements of the AMS matrix
  (\code{method=pairwise}), or alternatively, the columns similarity
  (\code{method=columnwise}), using K-means clustering (see
  function \sQuote{kmeans()}).

  Compared to the conventional cross-correlation matrix (see function
  \code{\link{dccm}}) the \dQuote{geostas} approach provide
  functionality to also detect domains involved in rotational
  motions (i.e. two  atoms located on opposite sides of a rotating
  domain will appear as anti-correlated in the cross-correlation matrix,
  but should obtain a high similarity coefficient in the AMS matrix). 

  See examples for more details.
}
\note{
  The current implementation in Bio3D uses a different fitting and
  clustering approach than the original Java implementation. The results
  will therefore differ.
}
\value{
  Returns a list object with the following components:
  \item{grps }{ a numeric vector containing the domain assignment per
    residue. }
  \item{amsm }{ a numeric matrix of atomic movement similarity (AMSM). }
  \item{fit.inds }{ a numeric vector of xyz indices used for fitting. }
}
\references{
  Romanowska, J. et al. (2012) \emph{JCTC} \bold{8}, 2588--2599.
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Julia Romanowska and Lars Skjaerven }
\seealso{
  \code{\link{read.ncdf}}, \code{\link{read.dcd}},
  \code{\link{read.pdb}}, \code{\link{mktrj}}.
}
\examples{
## Fetch stucture
pdb <- read.pdb("1hel")

## Calculate (vibrational) normal modes
modes <- nma(pdb)

## Mode trajectories
trj <- rbind(mktrj(modes, mode=7),
             mktrj(modes, mode=8),
             mktrj(modes, mode=9))

## Find domains
gs <- geostas(trj)

plot.dccm(gs$amsm, at=seq(0, 1, 0.1), main="AMSM with Domain Assignment ",
          margin.segments=gs$grps)


## Make a new PDB file with chain ID corresponding to domains
## write.pdb(pdb=NULL, xyz=modes$xyz, chain=gs$grps)

## gs <- geostas(trj, amsm=gs$amsm, k=5)

\dontrun{

#### Trajectory example
## Read inn DCD trajectory file, fit coordinates
dcdfile <- system.file("examples/hivp.dcd", package = "bio3d")
trj <- read.dcd(dcdfile)
xyz <- fit.xyz(trj[1,], trj)

## Find domains
gs <- geostas(xyz, k=3, fit=FALSE)

## Principal component analysis 
pc.md <- pca.xyz(xyz)

## Visualize PCs with colored domains (chain ID)
a <- mktrj(pc.md, pc=1, chain=gs$grps)



#### NMR-ensemble example
## Read a multi-model PDB file 
pdb <- read.pdb("1d1d", multi=TRUE)
ca.inds <- atom.select(pdb, 'calpha')

## Find domains and write PDB
gs  <- geostas(pdb)

## Fit all frames to the 'first' domain
domain.inds <- atom2xyz(which(gs$grps==1))

xyz <- fit.xyz(pdb$xyz.models[1, ca.inds$xyz],
               pdb$xyz.models[, ca.inds$xyz],
               fixed.inds = domain.inds,
               mobile.inds = domain.inds)

write.pdb(xyz=xyz, chain=gs$grps)



#### X-ray ensemble GroEL subunits
# Define the ensemble PDB-ids
ids <- c("1sx4_[A,B,H,I]", "1xck_[A-B]", "1sx3_[A-B]", "4ab3_[A-B]")

# Download and split PDBs by chain ID
raw.files <- get.pdb(ids, path = "raw_pdbs", gzip = TRUE)
files <- pdbsplit(raw.files, ids, path = "raw_pdbs/split_chain/")

# Align and superimpose coordinates
pdbs <- pdbaln(files)
gs <- geostas(pdbs, k=4, fit=TRUE)

# Principal component analysis 
gaps.pos <- gap.inspect(pdbs$xyz)
xyz <- fit.xyz(pdbs$xyz[1, gaps.pos$f.inds],
               pdbs$xyz[, gaps.pos$f.inds],
               fixed.inds=gs$fit.inds,
               mobile.inds=gs$fit.inds)
pc.xray <- pca.xyz(xyz)

## Visualize PCs with colored domains (chain ID)
a <- mktrj(pc.xray, pc=1, chain=gs$grps)

}
}
\keyword{ analysis }