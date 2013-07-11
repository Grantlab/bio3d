\name{dssp}
\alias{dssp}
\title{ Secondary Structure Analysis with DSSP }
\description{
  Secondary structure assignment according to the method of Kabsch and Sander.
}
\usage{
dssp(pdb, exepath = "", resno=TRUE)
}
\arguments{
  \item{pdb}{ a structure object of class \code{"pdb"}, obtained from
    \code{\link{read.pdb}}. }
  \item{exepath}{ path to the \sQuote{DSSP} program on your system
    (i.e. the directory where \sQuote{DSSP} is stored). }
  \item{resno}{ logical, if TRUE output is in terms of residue numbers
    rather than residue index (position in sequence). } 
}
\details{
  This function calls the \sQuote{DSSP} program to define secondary
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
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.

  \sQuote{DSSP} is the work of Kabsch and Sander:  
  Kabsch and Sander (1983) \emph{Biopolymers.}  \bold{12}, 2577--2637.

  For information on obtaining \sQuote{DSSP}, see:\cr 
  \url{http://swift.cmbi.ru.nl/gv/dssp/}.
}
\author{ Barry Grant }
\note{
  A system call is made to the \sQuote{DSSP} program, which must be
  installed on your system and in the search path for executables.
}
\seealso{ \code{\link{read.pdb}}, \code{\link{stride}},
  \code{\link{torsion.pdb}}, \code{\link{torsion.xyz}},
  \code{\link{plot.bio3d}}
}
\examples{
\dontrun{
# Read a PDB file
pdb <- read.pdb("1bg2")
sse <- dssp(pdb)

# Helix data
sse$helix

# Precent SSE content
sum(sse$helix$length)/sum(pdb$calpha) * 100
sum(sse$sheet$length)/sum(pdb$calpha) * 100
}
}
\keyword{ utilities }
