\name{stride}
\alias{stride}
\title{ Secondary Structure Analysis with STRIDE }
\description{
  Secondary structure assignment according to the method of Frishman and Argos.
}
\usage{
stride(pdb, exefile = "stride", resno=TRUE)
}
\arguments{
  \item{pdb}{ a structure object of class '"pdb"', obtained from
   \code{\link{read.pdb}}. }
  \item{exefile}{ file path to the \sQuote{STRIDE} program on your system  (i.e.
   how is \sQuote{STRIDE} invoked). }
  \item{resno}{ logical, if TRUE output is in terms of residue numbers
    rather than residue index (position in sequence). } 
}
\details{
  This function calls the \sQuote{STRIDE} program to define secondary
  structure and psi and phi torsion angles.
}
\value{
  Returns a list with the following components:
  \item{helix}{ \sQuote{start}, \sQuote{end}, \sQuote{length},
    \sQuote{chain} and \sQuote{type} of helix, where start and end are
    residue numbers or residue index positions depending on the value of
    \dQuote{resno} input argument. }
  \item{sheet}{ \sQuote{start}, \sQuote{end} and \sQuote{length} of E
    type sse, where start and end are residue numbers \dQuote{resno}. }
  \item{turn}{ \sQuote{start}, \sQuote{end} and \sQuote{length} of T
    type sse, where start and end are residue numbers \dQuote{resno}. }
  \item{phi}{ a numeric vector of phi angles. }
  \item{psi}{ a numeric vector of psi angles. }
  \item{acc}{ a numeric vector of solvent accessibility. } 
  \item{sse}{ a character vector of secondary structure type per residue. }
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.

  \sQuote{STRIDE} is the work of Frishman and Argos:
  Frishman and Argos (1995) \emph{Proteins.}  \bold{3}, 566--579.

  For information on obtaining the \sQuote{STRIDE} program, see:\cr
  \url{http://webclu.bio.wzw.tum.de/stride/},
  or copy it from an installation of VMD.
}
\author{ Barry Grant }
\note{
  A system call is made to the \sQuote{STRIDE} program which must be
  installed on your system and in the search path for executables.
}
\seealso{ \code{\link{read.pdb}}, \code{\link{dssp}},
  \code{\link{torsion.pdb}}, \code{\link{torsion.xyz}}
}
\examples{
\dontrun{
# Read a PDB file
pdb <- read.pdb( "4q21" )
sse <- stride(pdb)

# Helix data
sse$helix

# Precent SSE content
sum(sse$helix$length)/sum(pdb$calpha) * 100
sum(sse$sheet$length)/sum(pdb$calpha) * 100
}
}
\keyword{ utilities }
