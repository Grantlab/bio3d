\name{fit.xyz}
\alias{fit.xyz}
\alias{rot.lsq}
\title{ Coordinate Superposition }
\description{
  Coordinate superposition with the Kabsch algorithm.
}
\usage{
fit.xyz(fixed, mobile,
        fixed.inds  = NULL,
        mobile.inds = NULL,
        verbose=FALSE,
        pdb.path = "", pdbext = "",
        outpath = "fitlsq/", het=FALSE, full.pdbs=FALSE, ...)

rot.lsq(xx, yy,
        xfit = rep(TRUE, length(xx)), yfit = xfit,
        verbose = FALSE)
}
\arguments{
  \item{fixed }{ numeric vector of xyz coordinates.}
  \item{mobile}{ numeric vector, numeric matrix, or an object with an
    \code{xyz} component containing one or more coordinate sets. }
  \item{fixed.inds}{ a vector of indices that selects the elements of
    \code{fixed} upon which fitting should be based.}
  \item{mobile.inds}{  a vector of indices that selects the elements
    of \code{mobile} upon which fitting should be based.}
  \item{full.pdbs}{ logical, if TRUE \dQuote{full} coordinate files
    (i.e. all atoms) are written to the location specified by
    \code{outpath}. }
  \item{het}{ logical, if TRUE \sQuote{HETATM} records from full PDB
    files are written to output superposed PDB files.  Only required if
    \code{full.pdbs} is TRUE. } 
  \item{pdb.path}{ path to the \dQuote{full} input PDB files.  Only
        required if \code{full.pdbs} is TRUE. }
  \item{pdbext}{ the file name extension of the input PDB files. }
  \item{outpath}{ character string specifing the output directory when
    \code{full.pdbs} is TRUE. }
  \item{xx}{ numeric vector corresponding to the moving \sQuote{subject}
    coordinate set. }
  \item{yy}{ numeric vector corresponding to the fixed \sQuote{target}
    coordinate set. }
  \item{xfit}{ logical vector with the same length as \code{xx},  with
    TRUE elements corresponding to the subset of positions upon which
    fitting is to be performed. }
  \item{yfit}{ logical vector with the same length as \code{yy},  with
    TRUE elements corresponding to the subset of positions upon which
    fitting is to be performed. }
  \item{verbose}{ logical, if TRUE more details are printed. }
   \item{\dots}{ other parameters for \code{\link{read.pdb}}. }
  
}
\details{
  The function \code{fit.xyz} is a wrapper for the function
  \code{rot.lsq}, which performs the actual coordinate superposition.
  The function \code{rot.lsq} is an implementation of the Kabsch
  algorithm (Kabsch, 1978) and evaluates the optimal rotation matrix
  to minimize the RMSD between two structures.

  Since the Kabsch algorithm assumes that the number of points are the
  same in the two input structures, care should be taken to ensure that
  consistent atom sets are selected with \code{fixed.inds} and
  \code{mobile.inds}.

  Optionally, \dQuote{full} PDB file superposition and output can be
  accomplished by setting \cr \code{full.pdbs=TRUE}.  In that case, the
  input (\code{mobile}) passed to \code{fit.xyz} should be a list object
  obtained with the function \code{\link{read.fasta.pdb}}, since the
  components \code{id}, \code{resno} and \code{xyz} are required to
  establish correspondences.  See the examples below.
}
\value{
  Returns moved coordinates.
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.

  Kabsch \emph{Acta Cryst} (1978) \bold{A34}, 827--828.
}
\author{ Barry Grant with \code{rot.lsq} contributions from Leo Caves }
\seealso{ \code{\link{rmsd}}, \code{\link{read.pdb}},
  \code{\link{read.fasta.pdb}}, \code{\link{read.dcd}} }
\examples{
\dontrun{
##--- Read an alignment & Fit aligned structures
aln  <- read.fasta(system.file("examples/kif1a.fa",package="bio3d"))
pdb.path = system.file("examples",package="bio3d")
pdbs <- read.fasta.pdb(aln, pdb.path = pdb.path, pdbext = ".ent")
}

#-- OR Read previously saved kinesin data & Fit aligned structures
data(kinesin)
attach(kinesin, warn.conflicts=FALSE)

gaps <- gap.inspect(pdbs$xyz)

xyz <- fit.xyz( fixed  = pdbs$xyz[1,],
               mobile = pdbs$xyz,
               fixed.inds  = gaps$f.inds,
               mobile.inds = gaps$f.inds )

# Use indices from 'core.find()' for improved fitting 
#core <- core.find(pdbs) # see ?core.find

# Superpose again this time outputing PDBs
\dontrun{
xyz <- fit.xyz( fixed = pdbs$xyz[1,],
               mobile = pdbs,
               fixed.inds  = gaps$f.inds,
               mobile.inds = gaps$f.inds,
               pdb.path = system.file("examples/",package="bio3d"),
               pdbext = ".ent",
               outpath = "rough_fit/",
               full.pdbs = TRUE)

##--- Fit two PDBs
A <- read.pdb(system.file("examples/d1bg2__.ent",package="bio3d"))
A.ind <- atom.select(A, "///256:269///CA/")

B <- read.pdb(system.file("examples/d2kin.1.ent",package="bio3d"))
B.ind <- atom.select(B, "///257:270///CA/")

xyz <- fit.xyz(fixed=A$xyz, mobile=B$xyz,
               fixed.inds=A.ind$xyz,
               mobile.inds=B.ind$xyz)

# Write out moved PDB
C <- B; C$xyz = xyz
write.pdb(pdb=C, file = "moved.pdb")
}

}
\keyword{ utilities }

