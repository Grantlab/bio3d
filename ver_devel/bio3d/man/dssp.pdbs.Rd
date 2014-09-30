\name{dssp.pdbs}
\alias{dssp.pdbs}
\title{ Secondary Structure Analysis of Aligned PDB Structures with DSSP }
\description{
  Secondary structure assignment of aligned PDB structures -
  according to the method of Kabsch and Sander (DSSP) 
}
\usage{
dssp.pdbs(pdbs)
}
\arguments{
  \item{pdbs}{ a list object of class \code{"pdbs"} (obtained with
    \code{\link{pdbaln}} or \code{\link{read.fasta.pdb}}). }
}
\details{
  This is a wrapper function for easy SSE assignment of aligned PDB
  structures. 
}
\value{
  Returns a character matrix of aligned SSE assignment corresponding to
  each structure in the pdbs object (row  wise).
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Lars Skjaerven }
\seealso{ \code{\link{dssp}}, \code{\link{pdbaln}} }
\examples{
\dontrun{
## Fetch PDB files and split to chain A only PDB files
ids <- c("1a70_A", "1czp_A", "1frd_A", "1fxi_A", "1iue_A", "1pfd_A")
raw.files <- get.pdb(ids, path = "raw_pdbs")
files <- pdbsplit(raw.files, ids, path = "raw_pdbs/split_chain")

## Sequence Alignement
pdbs <- pdbaln(files)

## SSE assignment
sse <- dssp.pdbs(pdbs)
}
}
\keyword{ utilities }
