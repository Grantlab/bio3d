\name{split.pdb}
\alias{split.pdb}
\title{ Split a PDB File Into Separate Files, One For Each Chain. }
\description{
  Split a Protein Data Bank (PDB) coordinate file into new separate
  files with one file for each chain.
}
\usage{
split.pdb(pdb.files, path = "split_chain/")
}
\arguments{
  \item{pdb.files}{ a character vector of PDB file names. }
  \item{path}{ output path for chain-split files. }
}
\details{
  This function will produce single chain PDB files from multi-chain
  input files. 
}
\value{
  Called for its effect.
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
\seealso{ \code{\link{read.pdb}}, \code{\link{atom.select}}, \code{\link{write.pdb}}. }
\examples{
\dontrun{
pdbfiles <- list.files()
split.pdb(pdbfiles)
list.files("split_chain")
}
}
\keyword{ utilities }
