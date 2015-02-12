\name{geostas}
\alias{geostas}
\alias{geostas.default}
\alias{geostas.xyz}
\alias{geostas.nma}
\alias{geostas.enma}
\alias{geostas.pdb}
\alias{geostas.pdbs}
\alias{amsm.xyz}
\alias{print.geostas}
\title{ GeoStaS Domain Finder }
\description{
  Identifies geometrically stable domains in biomolecules
}
\usage{
geostas(\dots)

\method{geostas}{default}(\dots)

\method{geostas}{xyz}(xyz, amsm = NULL, k = 3, pairwise = TRUE,
      clustalg = "kmeans", fit = TRUE, ncore = NULL, verbose=TRUE, \dots)

\method{geostas}{nma}(nma, m.inds = 7:11, verbose=TRUE, \dots)

\method{geostas}{enma}(enma, pdbs = NULL, m.inds = 1:5, verbose=TRUE, \dots)

\method{geostas}{pdb}(pdb, inds = NULL, verbose=TRUE, \dots)

\method{geostas}{pdbs}(pdbs, verbose=TRUE, \dots)

amsm.xyz(xyz, ncore = NULL)

\method{print}{geostas}(x, \dots)
}
\arguments{
  \item{...}{ arguments passed to and from functions, such as
    \code{\link{kmeans}}, and \code{\link{hclust}} which are called
    internally in \code{geostas.xyz}. }
  \item{xyz}{ numeric matrix of xyz coordinates as obtained e.g. by
    \code{\link{read.ncdf}}, \code{\link{read.dcd}}, or
    \code{\link{mktrj}}. }
  \item{amsm}{ a numeric matrix as obtained by
    \code{\link{amsm.xyz}} (convenient e.g. for re-doing only the
    clustering analysis of the \sQuote{AMSM} matrix). }
  \item{k }{ an integer scalar or vector with the desired number of
    groups. }
  \item{pairwise }{ logical, if TRUE use pairwise clustering of the atomic
    movement similarity matrix (AMSM), else columnwise. }
  \item{clustalg}{ a character string specifing the clustering
    algorithm. Allowed values are \sQuote{kmeans} and \sQuote{hclust}. }
  \item{fit}{ logical, if TRUE coordinate superposition on identified
    core atoms is performed prior to the calculation of the AMS
    matrix. } 
  \item{ncore }{ number of CPU cores used to do the calculation.
    \code{ncore>1} requires package \sQuote{parallel} installed. }
  \item{verbose}{ logical, if TRUE details of the geostas calculations
    are printed to screen. }
  \item{nma}{ an \sQuote{nma} object as obtained from function
    \code{\link{nma}}. Function \code{\link{mktrj}} is used internally
    to generate a trajectory based on the normal modes. }
  \item{m.inds}{ the mode number(s) along which trajectory should be
    made (see function \code{\link{mktrj}}). }
  \item{enma}{ an \sQuote{enma} object as obtained from function
    \code{\link{nma.pdbs}}. Function \code{\link{mktrj}} is used internally
    to generate a trajectory based on the normal modes. }
  \item{pdbs}{ a \sQuote{pdbs} object as obtained from function
    \code{\link{pdbaln}} or \code{\link{read.fasta.pdb}}. }
  \item{pdb}{ a \sQuote{pdb} object as obtained from function
    \code{\link{read.pdb}}. }
  \item{inds}{ a \sQuote{select} object as obtained from function
    \code{\link{atom.select}} giving the atomic indices at which the
    calculation should be based. By default the function will attempt to
    locate C-alpha atoms using function \code{\link{atom.select}}. }
  \item{x}{ a \sQuote{geostas} object as obtained from function
    \code{\link{geostas}}. }
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
  (\code{pairwise=TRUE}), or alternatively, the columns similarity
  (\code{pairwise=FALSE}), using either K-means (\code{\link{kmeans}})
  or hierarchical (\code{\link{hclust}}) clustering.

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
  Returns a list object of type \sQuote{geostas} with the following components:
  \item{amsm }{ a numeric matrix of atomic movement similarity (AMSM). }
  \item{fit.inds }{ a numeric vector of xyz indices used for fitting. }
  \item{grps }{ a numeric vector containing the domain assignment per
    residue. }
  \item{atomgrps }{ a numeric vector containing the domain assignment per
    atom (only provided for \code{geostas.pdb}). }
  \item{inds }{ a list of atom \sQuote{select} objects with indices to
    corresponding to the identified domains. }
}
\references{
  Romanowska, J. et al. (2012) \emph{JCTC} \bold{8}, 2588--2599.
  Skjaerven, L. et al. (2014) \emph{BMC Bioinformatics} \bold{15}, 399.
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Julia Romanowska and Lars Skjaerven }
\seealso{
  \code{\link{plot.geostas}}, 
  \code{\link{read.pdb}}, \code{\link{mktrj}},
  \code{\link{read.ncdf}}, \code{\link{read.dcd}},
  \code{\link{nma}}, \code{\link{dccm}}.
}
\examples{
\donttest{

#### NMR-ensemble example
## Read a multi-model PDB file 
pdb <- read.pdb("1d1d", multi=TRUE)

## Find domains and write PDB
gs  <- geostas(pdb, fit=TRUE)

## Plot a atomic movement similarity matrix
plot.geostas(gs, contour=FALSE)

## Fit all frames to the 'first' domain
domain.inds <- gs$inds[[1]]

xyz <- pdbfit(pdb, inds=domain.inds)

#write.pdb(pdb, xyz=xyz, chain=gs$atomgrps)

}

\dontrun{
#### NMA example
## Fetch stucture
pdb <- read.pdb("1crn")

## Calculate (vibrational) normal modes
modes <- nma(pdb)

## Find domains
gs <- geostas(modes, k=2)

## Write NMA trajectory with domain assignment
mktrj(modes, mode=7, chain=gs$grps)

## Redo geostas domain clustering 
gs <- geostas(modes, amsm=gs$amsm, k=5)




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
mktrj(pc.md, pc=1, chain=gs$grps)




#### X-ray ensemble GroEL subunits
# Define the ensemble PDB-ids
ids <- c("1sx4_[A,B,H,I]", "1xck_[A-B]", "1sx3_[A-B]", "4ab3_[A-B]")

# Download and split PDBs by chain ID
raw.files <- get.pdb(ids, path = "raw_pdbs", gzip = TRUE)
files <- pdbsplit(raw.files, ids, path = "raw_pdbs/split_chain/")

# Align structures
pdbs <- pdbaln(files)

# Find domains
gs <- geostas(pdbs, k=4, fit=TRUE)

# Superimpose to core region
pdbs$xyz <- pdbfit(pdbs, inds=gs$fit.inds)

# Principal component analysis 
pc.xray <- pca(pdbs)

# Visualize PCs with colored domains (chain ID)
mktrj(pc.xray, pc=1, chain=gs$grps)


##- Same, but more manual approach 
gaps.pos <- gap.inspect(pdbs$xyz)

# Find core region
core <- core.find(pdbs)

# Fit to core region
xyz <- fit.xyz(pdbs$xyz[1, gaps.pos$f.inds],
               pdbs$xyz[, gaps.pos$f.inds],
               fixed.inds=core$xyz,
               mobile.inds=core$xyz)

# Find domains
gs <- geostas(xyz, k=4, fit=FALSE)

# Perform PCA
pc.xray <- pca.xyz(xyz)

# Make trajectory
mktrj(pc.xray, pc=1, chain=gs$grps)

}
}
\keyword{ analysis }
