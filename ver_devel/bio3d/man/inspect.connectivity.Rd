\name{inspect.connectivity}
\alias{inspect.connectivity}
\title{ Check the Connectivity of Protein Structures }
\description{
  Investigate protein coordinates to determine if the structure
  has missing residues. 
}
\usage{
inspect.connectivity(pdbs, cut=4.)
}
\arguments{
  \item{pdbs}{ an object of class \code{3daling} as obtained from
    function \code{pdbaln} or \code{read.fasta.pdb}; a xyz matrix
    containing the cartesian coordinates of C-alpha atoms; or a
    \sQuote{pdb} object as obtained from function \code{read.pdb}. }
  \item{cut }{ cutoff value to determine residue connectvitiy. }
}
\details{
  Utility function for checking if the PDB structures in a
  \sQuote{3dalign} object contains missing residues inside the
  structure. 

}
\value{
  Returns a vector. 
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Lars Skjaerven }
\seealso{ \code{\link{dm}}, \code{\link{gap.inspect}} }
\examples{
\dontrun{
## Fetch PDB files and split to chain A only PDB files
ids <- c("1a70_A", "1czp_A", "1frd_A", "1fxi_A", "1iue_A", "1pfd_A")
raw.files <- get.pdb(ids, path = "raw_pdbs")
files <- pdbsplit(raw.files, ids, path = "raw_pdbs/split_chain")

## Sequence Alignement, and connectivity check
pdbs <- pdbaln(files)

cons <- inspect.connectivity(pdbs)

## omit files with missing residues
files = files[cons]
}
}
\keyword{ analysis }