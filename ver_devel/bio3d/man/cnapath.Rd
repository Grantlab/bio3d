\name{cnapath}
\alias{cnapath}
\alias{summary.cnapath}
\alias{print.cnapath}
\title{ Suboptimal Path Analysis for Correlation Networks }
\description{
  Find k shortest paths between a pair of nodes, source and sink, in a correlation 
  network. 
}
\usage{
cnapath(cna, from, to=NULL, k=10, collapse=TRUE, ncore=NULL, \dots)
\method{summary}{cnapath}(object, \dots, pdb = NULL, label = NULL, col = NULL,
   plot = FALSE, concise = FALSE, cutoff = 0.1, normalize = TRUE, weight = FALSE)
\method{print}{cnapath}(x, \dots)
\method{plot}{cnapath}(x, \dots)
}
\arguments{
  \item{cna}{ A \sQuote{cna} object obtained from \code{\link{cna}}. }
  \item{from}{ Integer vector or matrix indicating node id(s) of source. If is matrix
       and \code{to} is NULL, the first column represents source and the second sink. }
  \item{to}{ Integer vector indicating node id(s) of sink. All combinations of 
       \code{from} and \code{to} values will be used as source/sink pairs. }
  \item{k}{ Integer, number of suboptimal paths to identify. }
  \item{collapse}{ Logical, if TRUE results from all source/sink pairs are merged with
       a single \sQuote{cnapath} object returned. }
  \item{ncore}{ Number of CPU cores used to do the calculation. 
       By default (NULL), use all detected CPU cores. }

  \item{object}{ A \sQuote{cnapath} class of object obtained from 
       \code{cnapath}. Multiple \sQuote{object} input is allowed for 
       comparing paths from different networks. }
  \item{pdb}{ A \sQuote{pdb} class of object obtained from \code{\link{read.pdb}} 
       and is used as the reference for node residue ids (in \code{summary.cnapath}) or 
       for molecular visulaization with VMD (in \code{vmd.cnapath}). }
  \item{label}{ Character, label for paths identified from different networks. }
  \item{col}{ Colors for plotting statistical results for paths identified 
       from different networks. }
  \item{plot}{ Logical, if TRUE path length distribution and node degeneracy will be plotted. }
  \item{concise}{ Logical, if TRUE only \sQuote{on path} residues will be displayed in the
       node degeneracy plot. }
  \item{cutoff}{ Numeric, nodes with node degeneracy larger than \code{cutoff} are
       shown in the output. }
  \item{normalize}{ Logical, if TRUE node degeneracy is divided by the total 
       (weighted) number of paths. }
  \item{weight}{ Logical, if TRUE each path is weighted by path length in calculating the
    node degeneracty. }
  \item{x}{ A 'cnapath' class object as obtained from function \code{cna}.
  }
  \item{\dots}{ Additional arguments passed to igraph function 
       \code{\link[igraph:distances]{get.shortest.paths}} (in 
       the function \code{cnapath}), passed to \code{summary.cnapath} 
       (in \code{print.cnapath}), as additional paths for 
       comparison (in \code{summary.cnapath}). }
}
\value{
  The function \code{cnapath} returns a \sQuote{cnapath} class of list containing 
  following three components:
  \item{path}{ a list object containing all identified suboptimal paths.  
      Each entry of the list is a sequence of node ids for the path. }
  \item{epath}{ a list object containing all identified suboptimal paths. 
      Each entry of the list is a sequence of edge ids for the path. }
  \item{dist}{ a numeric vector of all path lengths. }

  The function \code{summary.cnapath} returns a matrix of (normalized)
  node degeneracy for \sQuote{on path} residues. 
}
\references{
  Yen, J.Y. (1971) \emph{Management Science} \bold{17}, 712--716.
}
\author{ Xin-Qiu Yao }
\seealso{
  \code{\link{cna}}, \code{\link{cna.dccm}},
  \code{\link{vmd.cna}}, \code{\link{vmd.cnapath}}, 
  \code{\link[igraph:distances]{get.shortest.paths}}. }
\examples{
\donttest{
# Redundant testing excluded

if (!requireNamespace("igraph", quietly = TRUE)) {
   message('Need igraph installed to run this example')
} else {

attach(transducin)
inds = match(c("1TND_A", "1TAG_A"), pdbs$id)

npdbs <- trim(pdbs, row.inds=inds)
gaps.res <- gap.inspect(npdbs$ali)

modes <- nma(npdbs)
cij <- dccm(modes)
net <- cna(cij, cutoff.cij=0.3)

# get paths
pa1 <- cnapath(net[[1]], from = 314, to=172, k=50)
pa2 <- cnapath(net[[2]], from = 314, to=172, k=50)

# print the information of a path
pa1

# print two paths simultaneously
pas <- list(pa1, pa2)
names(pas) <- c("GTP", "GDP")
print.cnapath(pas)

# Or, for the same effect,
# summary(pa1, pa2, label=c("GTP", "GDP"))

# replace node numbers with residue name and residue number in the PDB file
pdb <- read.pdb("1tnd")
pdb <- trim.pdb(pdb, atom.select(pdb, chain="A", resno=npdbs$resno[1, gaps.res$f.inds]))
print.cnapath(pas, pdb=pdb)

# plot path length distribution and node degeneracy
print.cnapath(pas, pdb = pdb, col=c("red", "darkgreen"), plot=TRUE)

# View paths in 3D molecular graphic with VMD
#vmd.cnapath(pa1, pdb, launch = TRUE)
#vmd.cnapath(pa1, pdb, colors = 7, launch = TRUE)
#vmd.cnapath(pa1, pdb, spline=TRUE, colors=c("pink", "red"), launch = TRUE)
#pdb2 <- read.pdb("1tag")
#pdb2 <- trim.pdb(pdb2, atom.select(pdb2, chain="A", resno=npdbs$resno[2, gaps.res$f.inds]))
#vmd.cnapath(pa2, pdb2, launch = TRUE)

detach(transducin)

}
}
}
\keyword{ utilities }
