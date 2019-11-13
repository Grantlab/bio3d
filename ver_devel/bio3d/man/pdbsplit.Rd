\name{pdbsplit}
\alias{pdbsplit}
\title{ Split a PDB File Into Separate Files, One For Each Chain. }
\description{
  Split a Protein Data Bank (PDB) coordinate file into new separate
  files with one file for each chain.
}
\usage{
pdbsplit(pdb.files, ids = NULL, path = "split_chain", overwrite=TRUE,
         verbose = FALSE, mk4=FALSE, ncore = 1, progress = NULL, \dots)
}
\arguments{
  \item{pdb.files}{ a character vector of PDB file names. }
  \item{ids}{ a character vector of PDB and chain identifiers
    (of the form: \sQuote{pdbId_chainId}, e.g. \sQuote{1bg2_A}). 
    Used for filtering chain IDs for output (in the above example only 
    chain A would be produced). }
  \item{path}{ output path for chain-split files. }
  \item{overwrite}{ logical, if FALSE the PDB structures will not be
    read and written if split files already exist. }
  \item{verbose}{ logical, if TRUE details of the PDB header and
    chain selections are printed. }
  \item{mk4}{ logical, if TRUE output filenames will use only the first
    four characters of the input filename (see \code{basename.pdb} for details).  }   
  \item{ncore}{ number of CPU cores used for the calculation.
    \code{ncore>1} requires package \sQuote{parallel} be installed. }
  \item{progress}{ progress bar for use with shiny web app. }
  \item{...}{ additional arguments to \code{read.pdb}. Useful e.g. for
    parsing multi model PDB files, including ALT records etc. in the
    output files. }
}
\details{
  This function will produce single chain PDB files from multi-chain
  input files. By default all separate filenames are returned. To return 
  only a subset of select chains the optional input \sQuote{ids} can be 
  provided to filter the output (e.g. to fetch only chain C, of a PDB 
  object with additional chains A+B ignored). See examples section for 
  further details.

  Note that multi model atom records will only split into individual 
  PDB files if \code{multi=TRUE}, else they are omitted. See examples.
}
\value{
  Returns a character vector of chain-split file names.
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
\dontrun{
  ## Save separate PDB files for each chain of a local or on-line file
  pdbsplit( get.pdb("2KIN", URLonly=TRUE) )


  ## Split several PDBs by chain ID and multi-model records
  raw.files <- get.pdb( c("1YX5", "3NOB") , URLonly=TRUE)
  chain.files <- pdbsplit(raw.files,  path=tempdir(), multi=TRUE)
  basename(chain.files)


  ## Output only desired pdbID_chainID combinations
  ## for the last entry (1f9j), fetch all chains
  ids <- c("1YX5_A", "3NOB_B", "1F9J")
  raw.files <- get.pdb( ids , URLonly=TRUE)
  chain.files <- pdbsplit(raw.files, ids, path=tempdir())
  basename(chain.files)
}
}
\keyword{ utilities }
