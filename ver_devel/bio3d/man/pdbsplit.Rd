\name{pdbsplit}
\alias{pdbsplit}
\title{ Split a PDB File Into Separate Files, One For Each Chain. }
\description{
  Split a Protein Data Bank (PDB) coordinate file into new separate
  files with one file for each chain.
}
\usage{
pdbsplit(pdb.files, ids = NULL, path = "split_chain", multi=FALSE)
}
\arguments{
  \item{pdb.files}{ a character vector of PDB file names. }
  \item{ids}{ a character vector of PDB and chain identifiers
    (\sQuote{pdbId_chainId}). Used for filtering chain IDs for output. }
  \item{path}{ output path for chain-split files. }
  \item{multi}{ logical, if TRUE atom records for all models in
    multi-model files are splitted. }
}
\details{
  This function will produce single chain PDB files from multi-chain
  input files. By default all separate filenames are returned. To avoid
  this behaviour \sQuote{ids} can be provided to filter the
  output (e.g. to fetch only chain C, of a PDB object with additional
  chains A+B). See examples for details.

  Multi model atom records will be splitted into individual PDB files if
  \code{multi=TRUE}, else they are omitted.
}
\value{
  Returns a character vector of file names.
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
  
  For a description of PDB format (version3.3) see:\cr
  \url{http://www.wwpdb.org/documentation/format33/v3.3.html}.
}
\author{ Barry Grant }
\seealso{ \code{\link{read.pdb}}, \code{\link{atom.select}},
  \code{\link{write.pdb}}, \code{\link{get.pdb}}. }
\examples{
  ## Fetch PDBs
  pdb.ids <- c("1YX5", "3NOB", "1P3Q", "2D3G", "1F9J")
  outdir <- paste(tempdir(), "/raw_pdbs", sep="")
  raw.files <- get.pdb(pdb.ids, path = outdir)

  ## Split PDBs by chain ID and multi-model records
  all.files <- pdbsplit(raw.files,
            path =paste(outdir, "/split_chain", sep=""), multi=TRUE )
  
  list.files(paste(outdir, "/split_chain", sep=""))

\dontrun{
  ## desired pdbID-chainID combinations
  ids <- c("1YX5_A", "3NOB_B", "1P3Q_Q", "2D3G_A", "1F9J_A")

  ## Same as above, but return only selected pdbID-chainID combinations
  sel.files <- pdbsplit(raw.files, ids,
            path =paste(outdir, "/split_chain", sep=""), multi=TRUE )
}
}
\keyword{ utilities }
