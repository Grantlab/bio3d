\name{basename.pdb}
\alias{basename.pdb}
\title{
Manipulate PDB File Names
}
\description{
  Removes all of the path up to and including the last 
  path separator (if any) and the final \sQuote{.pdb} extension.
}
\usage{
basename.pdb(x, mk4 = FALSE, ext=".pdb")
}
\arguments{
  \item{x}{ character vector of PDB file names, containing path and extensions.}
  \item{mk4}{ logical, if TRUE the output will be truncated to the first 4 
    characters of the basename. This is frequently convenient for matching RCSB 
    PDB identifier conventions (see examples below). }
  \item{ext}{ character, specifying the file extension,
    e.g. \sQuote{.pdb} or \sQuote{.mol2}. }
}
\details{
This is a simple utility function for the common task of PDB file name manipulation. It is used internally in several bio3d functions and van be thought of as basename for PDB files. 
}
\value{
A character vector of the same length as the input \sQuote{x}.

Paths not containing any separators are taken to be in the current directory.

If an element of input is \sQuote{x} is \sQuote{NA}, so is the result.
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Barry Grant }
\seealso{
	\code{\link{basename}}, \code{\link{dirname}}
}
\examples{
	basename.pdb("/somedir/somewhere/1bg2_myfile.pdb")
	basename.pdb("/somedir/somewhere/1bg2_myfile.pdb", TRUE)
}
\keyword{ utilities }
