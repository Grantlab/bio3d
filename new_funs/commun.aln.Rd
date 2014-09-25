\name{commun.aln}
\alias{commun.aln}
\title{
  Compare communities in an ensemble of networks and re-color for equivalent communities
}
\description{
  This function facilitates the comparison between networks with same number 
  of nodes but different edges. It calculates the dissimilarity between all 
  communities and clusters communities with \sQuote{hclust} funciton. In each
  cluster, 0 or 1 community per network is included. Communities in the same
  cluster are then modified to have the same color.
}
\usage{
  commun.aln(nets, col = vmd.colors())
}
\arguments{
  \item{nets}{ a list containing an ensemble of \sQuote{cna} objects
    obtained from \sQuote{cna} function. }
  \item{col}{ a character vector containing colors for aligned community clusters. }
}
\value{
  Returns a list with the following components:
  \item{nets }{a list of \sQuote{cna} object after re-coloring. }
  \item{distmat }{dissimilarity matrix of bewteen all community pairs. }
  \item{hc }{returned value of calling \sQuote{hclust}. }
  \item{grps }{cluster numbers for all communities in all networks. }
}
\note{
  Only colors of nodes in the network and community.network are changed. 
  Membership numbering remains the original!
}
\details{ 
  Called for its effect.
}  
\seealso{
   \code{\link{cna}}, \code{\link{cna.cijs}}, \code{\link{dccm.enma}}, 
   \code{\link{nma.pdbs}}
}
\examples{

\dontrun{
data(transducin)
attach(transducin)

gaps.pos <- gap.inspect(pdbs$xyz)
modes <- nma.pdbs(pdbs, full=TRUE)
dccms <- dccm.enma(modes)
cijs <- dccms$all.dccm

nets <- cna.cijs(cijs)
nets.aln <- commun.aln(nets)$nets

# plot all-residue and coarse-grained (community) networks
pdb <- pdbs2pdb(pdbs, inds=1, rm.gaps=TRUE)[[1]]
op <- par(no.readonly=TRUE)
par(mar=c(0.1, 0.1, 0.1, 0.1))
layout(matrix(1:6, 2, 3, byrow=TRUE))

# before alignment
lapply(nets[1:3], function(x, pdb) {
    layout <- layout.cna(x, pdb, k=3)[,1:2]
    plot(x, layout=layout)}, pdb=pdb)

# after alignment
lapply(nets.aln[1:3], function(x, pdb) {
    layout <- layout.cna(x, pdb, k=3)[,1:2]
    plot(x, layout=layout)}, pdb=pdb)

layout(1)
par(op)
}
}
\keyword{analysis}
