\name{cna.cijs}
\alias{cna.cijs}
\title{
  Ensemble Protein Dynamic Correlation Network Construction and Community Analysis
}
\description{
  This function applies \sQuote{cna} function on an ensemble of correlation 
  matrix, as obtained from \sQuote{nma.pdbs}. The function builds both 
  residue-based and community-based undirected weighted
  network graphs for each matrix, and return a list of \sQuote{cna} objects.
  Community detection/clustering is performed on the initial residue based 
  network to determine the community organization and network structure
  of the community based network.
}
\usage{
  cna.cijs(cijs, \dots, ncore = NULL)
}
\arguments{
  \item{cijs}{ A numeric array with 3 dimensions (nXnXm) containing atomic
    correlation values, where "n" is the residue number and "m" the 
    system/frame number, or a \sQuote{list} object with each cell containing
    nXn matrix. The matrix elements should be in between 0 and 1 
    (atomic correlations). See
     \sQuote{dccm} function in bio3d package for further details. }
  \item{\dots}{ additional arguments to \code{cna}. }
  \item{ncore}{ number of CPU cores used to do the calculation. By default,
     all cores will be used}
}
\value{
  Returns a list of \sQuote{cna} object.
}
\details{ 
  Called for its effect.
}  
\seealso{
   \code{\link{cna}}, \code{\link{dccm.enma}}, \code{\link{nma.pdbs}}
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

# plot all-residue and coarse-grained (community) networks
pdb <- pdbs2pdb(pdbs, inds=1, rm.gaps=TRUE)[[1]]
op <- par(no.readonly=TRUE)
par(mar=c(0.1, 0.1, 0.1, 0.1))
layout(matrix(1:6, 2, 3, byrow=TRUE))

lapply(nets[1:3], function(x, pdb) {
   layout <- layout.cna(x, pdb, k=3, full=TRUE)[,1:2]
   plot(x$network, vertex.size=5, vertex.label=NA, 
   layout=layout)}, pdb=pdb)
lapply(nets[1:3], function(x, pdb) {
    layout <- layout.cna(x, pdb, k=3)[,1:2]
    plot(x, layout=layout)}, pdb=pdb)

layout(1)
par(op)
}
}
\keyword{analysis}
