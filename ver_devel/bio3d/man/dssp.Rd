\name{dssp}
\alias{dssp}
\alias{dssp.pdb}
\alias{dssp.pdbs}
\alias{dssp.xyz}
\alias{stride}
\alias{print.sse}
\title{ Secondary Structure Analysis with DSSP or STRIDE }
\description{
  Secondary structure assignment according to the method of Kabsch and Sander (DSSP)
  or the method of Frishman and Argos (STRIDE).
}
\usage{
dssp(\dots)

\method{dssp}{pdb}(pdb, exefile = "dssp", resno=TRUE, full=FALSE, verbose=FALSE, \dots)

\method{dssp}{pdbs}(pdbs, \dots)

\method{dssp}{xyz}(xyz, pdb, \dots)

stride(pdb, exefile = "stride", resno=TRUE)

\method{print}{sse}(x, \dots)
}
\arguments{
  \item{pdb}{ a structure object of class \code{"pdb"}, obtained from
    \code{\link{read.pdb}}. }
  \item{exefile}{ file path to the \sQuote{DSSP} or \sQuote{STRIDE} program 
    on your system (i.e. how is \sQuote{DSSP} or \sQuote{STRIDE} invoked). }
  \item{resno}{ logical, if TRUE output is in terms of residue numbers
    rather than residue index (position in sequence). }
  \item{full}{ logical, if TRUE bridge pairs and hbonds columns are
    parsed. }
  \item{verbose}{ logical, if TRUE \sQuote{DSSP} warning and error
    messages are printed. }
  \item{pdbs}{ a list object of class \code{"pdbs"} (obtained with
    \code{\link{pdbaln}} or \code{\link{read.fasta.pdb}}). }
  \item{xyz}{ a trajectory object of class \code{"xyz"}, obtained from
    \code{\link{read.ncdf}}, \code{\link{read.dcd}},
    \code{\link{read.crd}}. }
  \item{x}{ an \code{sse} object obtained from \code{\link{dssp.pdb}} or 
    \code{\link{stride}}. }
  \item{...}{ additional arguments to and from functions. }
  }
\details{
  This function calls the \sQuote{DSSP} or \sQuote{STRIDE} program to define 
  secondary structure and psi and phi torsion angles.
}
\value{
  Returns a list with the following components:
  \item{helix}{ \sQuote{start}, \sQuote{end}, \sQuote{length},
    \sQuote{chain} and \sQuote{type} of helix, where start and end are
    residue numbers or residue index positions depending on the value of
    \dQuote{resno} input argument. }
  \item{sheet}{ \sQuote{start}, \sQuote{end} and \sQuote{length} of E
    type sse, where start and end are residue numbers \dQuote{resno}. }
  \item{turn}{ \sQuote{start}, \sQuote{end} and \sQuote{length} of T
    type sse, where start and end are residue numbers \dQuote{resno}. }
  \item{phi}{ a numeric vector of phi angles. }
  \item{psi}{ a numeric vector of psi angles. }
  \item{acc}{ a numeric vector of solvent accessibility. }
  \item{sse}{ a character vector of secondary structure type per residue. }
  \item{hbonds}{ a 10 or 16 column matrix containing the bridge pair
    records as well as backbone NH-->O and O-->NH H-bond records. 
    (Only available for \code{\link{dssp}} }
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.

  \sQuote{DSSP} is the work of Kabsch and Sander:
  Kabsch and Sander (1983) \emph{Biopolymers.}  \bold{12}, 2577--2637.

  For information on obtaining \sQuote{DSSP}, see:\cr
  \url{http://swift.cmbi.ru.nl/gv/dssp/}.

  \sQuote{STRIDE} is the work of Frishman and Argos:
  Frishman and Argos (1995) \emph{Proteins.}  \bold{3}, 566--579.

  For information on obtaining the \sQuote{STRIDE} program, see:\cr
  \url{http://webclu.bio.wzw.tum.de/stride/},
  or copy it from an installation of VMD.
}
\author{ Barry Grant, Lars Skjaerven (dssp.pdbs) }
\note{
  A system call is made to the \sQuote{DSSP} or \sQuote{STRIDE} program, 
  which must be installed on your system and in the search path for executables.

  For the \code{hbonds} list component the column names can be
  used as a convenient means of data access, namely:\cr
  Bridge pair 1            \dQuote{BP1},\cr
  Bridge pair 2            \dQuote{BP2},\cr
  Backbone H-bond (NH-->O) \dQuote{NH-O.1},\cr
  H-bond energy of NH-->O  \dQuote{E1},\cr
  Backbone H-bond (O-->NH) \dQuote{O-HN.1},\cr
  H-bond energy of O-->NH  \dQuote{E2},\cr
  Backbone H-bond (NH-->O) \dQuote{NH-O.2},\cr
  H-bond energy of NH-->O  \dQuote{E3},\cr
  Backbone H-bond (O-->NH) \dQuote{O-HN.2},\cr
  H-bond energy of O-->NH  \dQuote{E4}.\cr

  If \sQuote{resno=TRUE} the following additional columns are included:\cr
  Chain ID of resno \dQuote{BP1}:    \dQuote{ChainBP1},\cr
  Chain ID of resno \dQuote{BP2}:    \dQuote{ChainBP2},\cr
  Chain ID of resno \dQuote{O-HN.1}: \dQuote{Chain1},\cr
  Chain ID of resno \dQuote{NH-O.2}: \dQuote{Chain2},\cr
  Chain ID of resno \dQuote{O-HN.1}: \dQuote{Chain3},\cr
  Chain ID of resno \dQuote{NH-O.2}: \dQuote{Chain4}.
}
\seealso{
  \code{\link{read.pdb}},
  \code{\link{torsion.pdb}}, \code{\link{torsion.xyz}},
  \code{\link{plot.bio3d}},
  
  \code{\link{read.ncdf}}, \code{\link{read.dcd}},
  \code{\link{read.prmtop}}, \code{\link{read.crd}},
}
\examples{
\dontrun{
##- PDB example
# Read a PDB file
pdb <- read.pdb("1bg2")
sse <- dssp(pdb)
sse2 <- stride(pdb)

## Short summary
sse
sse2

# Helix data
sse$helix

# Precent SSE content
sum(sse$helix$length)/sum(pdb$calpha) * 100
sum(sse$sheet$length)/sum(pdb$calpha) * 100



##- PDBs example
## Fetch PDB files and split to chain A only PDB files
ids <- c("1a70_A", "1czp_A", "1frd_A", "1fxi_A", "1iue_A", "1pfd_A")
raw.files <- get.pdb(ids, path = "raw_pdbs")
files <- pdbsplit(raw.files, ids, path = "raw_pdbs/split_chain")

## Sequence Alignement
pdbs <- pdbaln(files)

## SSE assignment
sse <- dssp.pdbs(pdbs)



##- XYZ Trajectory
pdb <- read.pdb("2mda", multi=TRUE)
dssp.xyz(pdb$xyz, pdb)

}
}
\keyword{ utilities }
