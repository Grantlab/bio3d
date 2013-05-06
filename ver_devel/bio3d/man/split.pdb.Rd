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
  
  For a description of PDB format (version3.3) see:\cr
  \url{http://www.wwpdb.org/documentation/format33/v3.3.html}.
}
\author{ Barry Grant }
\note{
  To Be Improved.
}
\seealso{ \code{\link{read.pdb}}, \code{\link{atom.select}}, \code{\link{write.pdb}}. }
\examples{
\dontrun{
#pdbfiles <- list.files()
pdbfiles <- get.pdb(c("1bg2", "1cz7"), path=tempdir())
split.pdb(pdbfiles, path=paste(tempdir(), "/split_chain/", sep=""))
list.files(paste(tempdir(),"/split_chain",sep=""))
}
}
\keyword{ utilities }
