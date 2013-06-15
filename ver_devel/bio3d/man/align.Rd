\name{align}
\alias{align}
\title{ 3D alignment }
\description{
  Performs a sequence and structural alignment of two PDB entities.
}
\usage{
align(fixed, mobile, fixed.inds=NULL, mobile.inds=NULL, write.pdbs=TRUE,
      max.cycles=10, cutoff=0.5, seqaln.file="aln.fa")
}
\arguments{
  \item{fixed}{ an object of class \code{pdb} as obtained from
    function \code{read.pdb}. }
  \item{mobile}{ an object of class \code{pdb} as obtained from
    function \code{read.pdb}. }
  \item{fixed.inds}{ atom and xyz coordinate indices obtained from \code{atom.select}
    that selects the elements of \code{fixed} upon which the calculation
    should be based.}
  \item{mobile.inds}{ atom and xyz coordinate indices obtained from \code{atom.select}
    that selects the elements of \code{mobile} upon which the calculation
    should be based.}
  \item{write.pdbs}{ logical, if TRUE the aligned structures are written
    to file.}
  \item{max.cycles}{ maximum number of refinement cycles. }
  \item{cutoff}{ standard deviation at which the fitting refinement stops.}
  \item{seqaln.file}{ name of \sQuote{FASTA} output file to which alignment should be
    written.}
}
\details{
  This function performs a sequence alignment followed by a structural
  alignment of the two PDB entities. Cycles of refinement steps of the
  structural alignment are performed to improve the fit by removing
  atoms with a high structural deviation. The primary purpose of the function
  is to allow rapid structural alignment for distant protein structures.

  The function reports the residues of \code{fixed} and \code{mobile}
  included in the final structural alignment, as well as the related RMSD values.
  
  This is largely a wrapper function calling the underlying function
  \code{seqaln}, \code{rot.lsq}, and \code{rmsd}.
}
\value{
  Returns a list with the following components:
  \item{a.inds}{ atom and xyz indices of \code{fixed}. }
  \item{b.inds}{ atom and xyz indices of \code{mobile}. }
  \item{rmsd}{ a numeric vector of RMSD values after cycles of
    refinement steps. }
}
\references{
   Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Lars Skjærven }
\seealso{ \code{\link{rmsd}}, \code{\link{rot.lsq}}, \code{\link{seqaln}} }
\examples{

     ## Stucture of PKA:
     a <- read.pdb("1cmk")

     ## Stucture of PKB:
     b <- read.pdb("2jdo")

     ## Align and fit b on to a:
     aln <- align(a, b)

     ## Should be the same as aln$rmsd
     rmsd(a$xyz, b$xyz, aln$a.inds$xyz, aln$b.inds$xyz, fit=TRUE)


     ## Align two subunits of GroEL (open and closed states)
     a <- read.pdb("1sx4", maxlines=70000)
     b <- read.pdb("1xck", maxlines=70000)

     ## Select chain A only
     a.inds <- atom.select(a, chain="A")
     b.inds <- atom.select(b, chain="A")

     ## Align and fit:
     aln <- align(a,b, a.inds, b.inds)

     ## Fit the entire complex based on returned indices
     xyz <- fit.xyz(a$xyz, b$xyz, aln$a.inds$xyz, aln$b.inds$xyz)
     write.pdb(b, xyz=xyz)

 
}
\keyword{ utilities }