\name{dssp.trj}
\alias{dssp.trj}
\title{ Secondary Structure Analysis of Trajectories with DSSP  }
\description{
  Secondary structure assignment according to the method of Kabsch and Sander.
}
\usage{
dssp.trj(pdb, trj, skip=1000, threshold=3, file.head="")
}
\arguments{
  \item{pdb}{ a structure object of class \code{"pdb"}, obtained from
    \code{\link{read.pdb}}. }
  \item{trj}{ a trajectory object of class \code{"trj"}, obtained from \code{\link{read.ncdf}},  \code{\link{read.dcd}}, \code{\link{read.crd}}. }
  \item{skip}{ a number indicating the frame frequency in the trajectory indicated; default = 1000. }
  \item{threshold}{ a number indicating the threshold to be used for the analysis: higher numbers will decrease the threshold, allowing to accept structures with smaller secondary structure differences, with respect of the reference \code{pdb} structure; default = 3}
  \item{file.head}{ a path where to save the structures with possible unfolding events. }
}
\details{
  This function calls the \sQuote{DSSP} program and calculates the secondary structure difference between the reference \code{pdb} and frames from the \code{trj} file indicated, according with the \code{skip} value indicated.
}
\value{ This function gives, as output, a comparative secondary structure analysis between different structures: particularly, the frame number under analysis will be printed on screen if the structure shows possible unfolding events and the frame structure will be saved, according with the \code{file.head} indicated. 
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
\seealso{ \code{\link{read.pdb}}, \code{\link{read.ncdf}},
  \code{\link{read.dcd}}, \code{\link{read.crd}},
  \code{\link{plot.bio3d}}, \code{\link{dssp}}
}
\examples{
\dontrun{
# Read a PDB file
pdb <- read.pdb(system.file("examples/hivp.pdb", package="bio3d"))
trj <- read.dcd(system.file("examples/hivp.dcd", package="bio3d"))

dssp.trj(pdb, trj, skip = 100, threshold = 1.5)

}
}
\keyword{ utilities }
