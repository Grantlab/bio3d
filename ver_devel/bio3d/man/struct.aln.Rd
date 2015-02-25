\name{struct.aln}
\alias{struct.aln}
\title{ Structure Alignment Of Two PDB Files }
\description{
  Performs a sequence and structural alignment of two PDB entities.
}
\usage{
struct.aln(fixed, mobile, fixed.inds=NULL, mobile.inds=NULL, 
           write.pdbs=TRUE, outpath = "fitlsq", prefix=c("fixed",
           "mobile"), max.cycles=10, cutoff=0.5, ... )
}
\arguments{
  \item{fixed}{ an object of class \code{pdb} as obtained from
    function \code{read.pdb}. }
  \item{mobile}{ an object of class \code{pdb} as obtained from
    function \code{read.pdb}. }
  \item{fixed.inds}{ atom and xyz coordinate indices obtained from
    \code{atom.select} that selects the elements of \code{fixed} upon
    which the calculation should be based.}
  \item{mobile.inds}{ atom and xyz coordinate indices obtained from
    \code{atom.select} that selects the elements of \code{mobile} upon
    which the calculation should be based.}
  \item{write.pdbs}{ logical, if TRUE the aligned structures are written
    to PDB files.}
  \item{outpath}{ character string specifing the output directory when
    \code{write.pdbs} is TRUE. }
  \item{prefix}{ a character vector of length 2 containing the filename
    prefix in which the fitted structures should be written. }
  \item{max.cycles}{ maximum number of refinement cycles.}
  \item{cutoff}{ standard deviation of the pairwise distances for
    aligned residues at which the fitting refinement stops.}
  \item{\dots}{ extra arguments passed to \code{seqaln} function.}
}
\details{
  This function performs a sequence alignment followed by a structural
  alignment of the two PDB entities. Cycles of refinement steps of the
  structural alignment are performed to improve the fit by removing
  atoms with a high structural deviation. The primary purpose of the
  function is to allow rapid structural alignment (and RMSD analysis)
  for protein structures with unequal, but related sequences.

  The function reports the residues of \code{fixed} and \code{mobile}
  included in the final structural alignment, as well as the related
  RMSD values.
  
  This function makes use of the underlying functions \code{seqaln},
  \code{rot.lsq}, and \code{rmsd}.
}
\value{
  Returns a list with the following components:
  \item{a.inds}{ atom and xyz indices of \code{fixed}. }
  \item{b.inds}{ atom and xyz indices of \code{mobile}. }
  \item{xyz}{ fitted xyz coordinates of \code{mobile}. }
  \item{rmsd}{ a numeric vector of RMSD values after each cycle of
    refinement. }
}
\references{
   Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Lars Skjarven }
\seealso{ \code{\link{rmsd}}, \code{\link{rot.lsq}},
  \code{\link{seqaln}}, \code{\link{pdbaln}} }
\examples{
\donttest{
# Needs MUSCLE installed - testing excluded

if(check.utility("muscle")) {

     ## Stucture of PKA:
     a <- read.pdb("1cmk")

     ## Stucture of PKB:
     b <- read.pdb("2jdo")

     ## Align and fit b on to a:
     path = file.path(tempdir(), "struct.aln")
     aln <- struct.aln(a, b, outpath = path, outfile = tempfile())
     
     ## Should be the same as aln$rmsd (when using aln$a.inds and aln$b.inds)
     rmsd(a$xyz, b$xyz, aln$a.inds$xyz, aln$b.inds$xyz, fit=TRUE)
     
     invisible( cat("\nSee the output files:", list.files(path, full.names = TRUE), sep="\n") )
}
}

\dontrun{
     ## Align two subunits of GroEL (open and closed states)
     a <- read.pdb("1sx4")
     b <- read.pdb("1xck")

     ## Select chain A only
     a.inds <- atom.select(a, chain="A")
     b.inds <- atom.select(b, chain="A")

     ## Align and fit:
     aln <- struct.aln(a,b, a.inds, b.inds)
}
}
\keyword{ utilities }
