\name{pdbs2pdb}
\alias{pdbs2pdb}
\title{ PDBs to PDB Converter }
\description{
  Convert a list of PDBs from an \code{"3dalign"} object to a list of
  \code{pdb} objects. 
}
\usage{
pdbs2pdb(pdbs, inds = NULL, rm.gaps = FALSE)
}
\arguments{
  \item{pdbs}{ a list of class \code{"3dalign"} containing PDB file
    data, as obtained from \code{read.fasta.pdb} or
    \code{pdbaln}. }
  \item{inds}{ a vector of indices that selects the PDB structures to
    convert. }
  \item{rm.gaps}{ logical, if TRUE atoms in gap containing columns are
    removed in the output \code{pdb} objects. }
}
\details{
  This function will generate a list of \code{pdb} objects from a
  \code{"3dalign"} class.

  See examples for more details/
}
\value{
  Returns a list of \code{pdb} objects.
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Lars Skjaerven }
\seealso{ \code{\link{read.pdb}}, \code{\link{pdbaln}},
  \code{\link{read.fasta.pdb}}. }
\examples{
\dontrun{
## Fetch PDBs
pdb.ids <- c("1YX5_B", "3NOB", "1P3Q_U")
#outdir <- paste(tempdir(), "/raw_pdbs", sep="")
outdir = "raw_pdbs"
raw.files <- get.pdb(pdb.ids, path = outdir)

## Split PDBs by chain ID and multi-model records
all.files <- pdbsplit(raw.files, pdb.ids, 
                path =paste(outdir, "/split_chain", sep=""))

## Align and fit
pdbs     <- pdbaln(all.files, fit=TRUE)

## Convert back to PDB objects
all.pdbs <- pdbs2pdb(pdbs)

## Access the first PDB object
## all.pdbs[[1]]

## Return PDB objects consisting of only
## atoms in non-gap positions
all.pdbs <- pdbs2pdb(pdbs, rm.gaps=TRUE)

}
}
\keyword{ utilities }
