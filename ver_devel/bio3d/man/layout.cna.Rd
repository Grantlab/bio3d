\name{layout.cna}
\alias{layout.cna}
\title{ Protein Structure Network Layout }
\description{
  Determine protein structure network layout in 2D and 3D from the geometric 
  center of each community.
}
\usage{
  layout.cna(x, pdb, renumber=TRUE, k=2, full=FALSE)
}
\arguments{
  \item{x}{ A protein structure network object as obtained from the \sQuote{cna} 
    function. }
  \item{pdb}{ A pdb class object as obtained from the \sQuote{read.pdb} 
    function. }
  \item{renumber}{ Logical, if TRUE the input \sQuote{pdb} will be re-numbered 
    starting at residue number one before community coordinate averages are 
    calculated. } 
  \item{k}{ A single element numeric vector between 1 and 3 specifying the 
    returned coordinate dimensions. }
  \item{full}{ Logical, if TRUE the full all-Calpha atom network coordinates will 
    be returned rather than the default clustered network community coordinates. }       

}
\details{
  This function calculates the geometric center for each community from
  the atomic position of it's Calpha atoms taken from a corresponding PDB 
  file. Care needs to be taken to ensure the PDB residue numbers and the
  community vector names/length match.

  The community residue membership are typically taken from the input network 
  object but can be supplied as a list object with 'x$communities$membership'.
}
\value{
 A numeric matrix of Nxk, where N is the number of communities and k the number of dimensions requested.
}
\author{ Guido Scarabelli and Barry Grant }

\seealso{ \code{\link{plot.cna}}, \code{\link[igraph:plot.communities]{plot.communities}},
  \code{\link[igraph:igraph.plotting]{igraph.plotting}}, 
  \code{\link[igraph:plot.igraph]{plot.igraph}}}
\examples{

# Load the correlation network
attach(hivp)

# Read the starting PDB file to determine atom correspondence
pdbfile <- system.file("examples/hivp.pdb", package="bio3d")
pdb <- read.pdb(pdbfile)

# Plot will be slow
xy <- plot.cna(net)
#plot3d.cna(net, pdb)

layout.cna(net, pdb, k=3)
layout.cna(net, pdb)

# can be used as input to plot.cna and plot3d.cna....
# plot.cna( net, layout=layout.cna(net, pdb) )
# plot3d.cna(net, pdb, layout=layout.cna(net, pdb, k=3))

}
\keyword{ utility }
