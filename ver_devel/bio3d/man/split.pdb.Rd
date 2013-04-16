\name{split.pdb}
\alias{split.pdb}
\title{ Split a PDB File Into Separate Files, One For Each Chain. }
\description{
  Split a Protein Data Bank (PDB) coordinate file into new separate
  files with one file for each chain.
}
\usage{
split.pdb(pdb.files, path = "split_chain/", multi=FALSE)
}
\arguments{
   \item{pdb.files}{ a character vector of PDB file names. }
  \item{path}{ output path for chain-split files. }
  \item{multi}{ logical, if TRUE atom records for all models in
    multi-model files are splitted. }
}
\details{
  This function will produce single chain PDB files from multi-chain
  input files.

  Multi model atom records will be splitted into individual PDB files if
  \code{multi=TRUE}, else they are omitted.
}
\value{
  Returns a character vector of file names.
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
  
  For a description of PDB format (version2.2) see:\cr
  \url{http://www.rcsb.org/pdb/file_formats/pdb/pdbguide2.2/guide2.2_frame.html}
.
}
\author{ Barry Grant }
\note{
  To Be Improved.
}
\seealso{ \code{\link{read.pdb}}, \code{\link{atom.select}},
  \code{\link{write.pdb}}, \code{\link{get.pdb}}. }
\examples{
  ## Fetch PDBs
  pdb.ids <- c("1YX5", "3NOB", "1P3Q", "2D3G","1F9J")
  raw.files <- get.pdb(pdb.ids, path = "raw_pdbs")

  ## Split PDBs by chain ID and multi-model records
  all.files <- split.pdb(raw.files, path ="raw_pdbs/split_chain", multi=TRUE )

  \dontrun{
  pdbfiles <- list.files()
  split.pdb(pdbfiles)
  list.files("split_chain")
  }
}
\keyword{ utilities }
