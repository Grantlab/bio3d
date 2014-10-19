\name{prune.cna}
\alias{prune.cna}
\title{ Prune A cna Network Object }
\description{
Remove nodes and their associated edges from a cna network graph.
}
\usage{
prune.cna(x, edges.min = 1, size.min = 1)
}
\arguments{
  \item{x}{ A protein network graph object as obtained from the
    \sQuote{cna} function.}
  \item{edges.min}{ A single element numeric vector specifying the 
     minimum number of edges that retained nodes should have. Nodes 
     with less than \sQuote{edges.min} will be pruned. }
  \item{size.min}{ A single element numeric vector specifying the 
     minimum node size that retained nodes should have. Nodes 
     with less composite residues than \sQuote{size.min} will be 
     pruned. }
}
\details{
 This function is useful for cleaning up cna network plots by removing, 
for example, small isolated nodes. The output is a new cna object minus 
the pruned nodes and their associated edges. Node naming is preserved.
}
\value{
A cna class object, see function \code{\link{cna}} for details.
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{
Barry Grant
}
\note{
Some improvements to this function are required, including a better effort 
to preserve the original community structure rather than calculating a 
new one. Also may consider removing nodes form the raw.network object 
that is returned also.
}
\seealso{
   \code{\link{cna}},  \code{\link{summary.cna}},
   \code{\link{view.cna}}, \code{\link{plot.cna}}
}
\examples{

# Load the correlation network
attach(hivp)

# Read the starting PDB file to determine atom correspondence
pdbfile <- system.file("examples/hivp.pdb", package="bio3d")
pdb <- read.pdb(pdbfile)

# Plot coarse grain network based on dynamically coupled communities
par(mfcol=c(1,2), mar=c(0,0,0,0))
plot.cna(net)

# Prune network
dnet <- prune.cna(net, edges.min = 1)
plot(dnet)
}
\keyword{ utility }
