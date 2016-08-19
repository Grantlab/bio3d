\name{pymol}
\alias{pymol}
\alias{pymol.pdbs}
\alias{pymol.nma}
\alias{pymol.pca}
\alias{pymol.modes}
\alias{pymol.dccm}
\title{ Biomolecular Visualization with PyMOL }
\description{
  Visualize Bio3D structure objects in PyMOL
}
\usage{
pymol(\dots)

\method{pymol}{pdbs}(pdbs, col=NULL, as="ribbon", file=NULL, type="script",
                     exefile="pymol", \dots)

\method{pymol}{nma}(\dots)

\method{pymol}{pca}(\dots)

\method{pymol}{modes}(modes, mode=NULL, file=NULL, scale=5, dual=FALSE,
                     type="script", exefile="pymol", \dots)

\method{pymol}{dccm}(dccm, pdb, file=NULL,
                     step=0.2, omit=0.2, radius = 0.15,
                     type="script", exefile="pymol", \dots)


}
\arguments{
  \item{pdbs}{ aligned C-alpha Cartesian coordinates as obtained with 
    \code{\link{read.fasta.pdb}} or \code{\link{pdbaln}}. }
  \item{col}{ a single element character vector specifying the coloring
    of the structures. Options are: \sQuote{index}, \sQuote{index2},
    \sQuote{gaps}, \sQuote{rmsf}.

    Special cases: Provide a \sQuote{core} object as obtained by
    \code{core.find} to color on the invariant core. Alternatively,
    provide a vector containing the color code for each structure in the
    \sQuote{pdbs} object. }
  \item{as}{ show as \sQuote{ribbon}, \sQuote{cartoon}, \sQuote{lines}, \sQuote{putty}. }
  \item{file}{ a single element character vector specifying the file name
    of the PyMOL session/script file. }
  \item{type}{ a single element character vector specifying the output
    type: \sQuote{script} generates a .pml script; \sQuote{session}
    generates a .pse session file; \sQuote{launch} launches pymol. }
  \item{exefile}{ file path to the \sQuote{PYMOL} program on your system (i.e.
    how is \sQuote{PYMOL} invoked). }
  \item{modes}{ an object of class \code{nma} or \code{pca} as obtained
    from functions \code{nma} or \code{pca.xyz}. }
  \item{mode}{ the mode number for which the vector field should be
    made. } 
  \item{scale}{ global scaling factor. }
  \item{dual}{ logical, if TRUE mode vectors are also drawn in both
    direction. }
  \item{dccm}{ an object of class \code{dccm} as obtained from
    function \code{\link{dccm}}. }
  \item{pdb}{ an object of class \code{pdb} as obtained from
    function \code{read.pdb} or a numerical vector of Cartesian
    coordinates. }
  \item{step}{ binning interval of cross-correlation coefficents. }
  \item{omit}{ correlation coefficents with values (0-omit, 0+omit) will
    be omitted from visualization. }
  \item{radius}{ numeric, radius of visualized correlation cylinders in
    PyMol. Alternatively, a matrix with the same dimesions as
    \code{dccm} can be provided, e.g. to draw cylinders with radii
    associated to the pairwise correlation value. }
  \item{\dots}{ arguments passed to function \code{pymol.modes} for
    \sQuote{nma} and \sQuote{pca} objects. }
}
\details{
  These functions provides a convenient approach for the visualization of
  Bio3D objects in PyMOL. See examples for more details.

  DCCM PyMOL visualization:
  This function generates a PyMOL (python) script that will draw colored
  lines between (anti)correlated residues. The PyMOL script file is
  stored in the working directory with filename \dQuote{R.py}.
  PyMOL will only be launched (and opened)  when using argument
  \sQuote{type='launch'}. Alternatively a PDB file with CONECT records
  will be generated (when argument \code{type='pdb'}). 
   
  For the PyMOL version, PyMOL CGO objects are generated - each object
  representing a range of correlation values (corresponding to the
  actual correlation values as found in the correlation
  matrix). E.g. the PyMOL object with name \dQuote{cor_-1_-08} would
  display all pairs of correlations with values between -1 and -0.8.

  NMA / PCA PyMOL vector field visualization:
  This function generates a PyMOL (python) script for drawing
  mode vectors on a PDB structure. The PyMOL script file is stored in
  the working directory with filename \dQuote{R.py}. 

}
\value{
  Called for its action
}
\references{
  Skjaerven, L. et al. (2014) \emph{BMC Bioinformatics} \bold{15}, 399.
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Lars Skjaerven }
\seealso{
  \code{view}
}
\examples{
\dontrun{

##- pymol with a 'pdbs' object
attach(transducin)

# build a pymol session containing all structures in the PDBs object
pymol(pdbs)

# color by invariant core (
# core <- core.find(pdbs)
pymol(pdbs, col=core)

# color by RMSF
pymol(pdbs, col="rmsf")

# color by clustering
rd <- rmsd(pdbs$xyz)
hc <- hclust(as.dist(rd))
grps <- cutree(hc, k=3)
pymol(pdbs, col=grps)



##- pymol with a 'dccm' object
## Fetch stucture
pdb <- read.pdb( system.file("examples/1hel.pdb", package="bio3d") )

## Calculate normal modes
modes <- nma(pdb)

## Calculate correlation matrix
cm <- dccm.nma(modes)

pymol(cm, modes$xyz)



##- pymol with a 'nma' or 'pca' object
pymol(modes, mode=7)


}
}
