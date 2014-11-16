\name{sse.bridges}
\alias{sse.bridges}
\title{ SSE Backbone Hydrogen Bonding }
\description{
  Determine backbone C=O to N-H hydrogen bonding in secondary
  structure elements.
}
\usage{
sse.bridges(sse, type="helix", hbond=TRUE, energy.cut=-1.0)
}
\arguments{
  \item{sse}{ an sse object as obtained with \code{dssp}. }
  \item{type}{ character string specifying \sQuote{helix} or \sQuote{sheet}. }
  \item{hbond}{ use hbond records in the dssp output. }
  \item{energy.cut}{ cutoff for the dssp hbond energy. }
}
\details{
  Simple functionality to parse the \sQuote{BP} and \sQuote{hbond}
  records of the DSSP output.

  Requires input from function \code{dssp} with arguments
  \code{resno=FALSE} and \code{full=TRUE}. 
}
\value{
  Returns a numeric matrix of two columns containing the residue ids
  of the paired residues.
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Lars Skjaerven }
\seealso{ \code{\link{read.pdb}}, \code{\link{dssp}} }
\examples{
\dontrun{
# Read a PDB file
pdb <- read.pdb( system.file("examples/1hel.pdb", package="bio3d") )
sse <- dssp(pdb, resno=FALSE, full=TRUE)

sse.bridges(sse, type="helix")
}
}
\keyword{ utilities }
