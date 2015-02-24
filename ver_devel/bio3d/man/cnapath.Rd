\name{cnapath}
\alias{cnapath}
\alias{summary.cnapath}
\alias{print.cnapath}
\alias{view.cnapath}
\title{ Suboptimal Path Analysis for Correlation Networks }
\description{
  Find k shortest paths between a pair of nodes, source and sink, in a correlation 
  network. 
}
\usage{
cnapath(cna, from, to, k = 10, ncore = NULL, \dots)
\method{summary}{cnapath}(object, \dots, pdb = NULL, label = NULL, col = NULL,
   plot = FALSE, concise = FALSE, cutoff = 0.1, normalize = TRUE)
\method{print}{cnapath}(x, \dots)
view.cnapath(x, pdb, out.prefix = "view.cnapath", spline = FALSE, 
   colors = c("blue", "red"), launch = FALSE, \dots)
}
\arguments{
  \item{cna}{ A \sQuote{cna} object obtained from \code{\link{cna}}. }
  \item{from}{ Integer, node id for the source. }
  \item{to}{ Integer, node id for the sink. }
  \item{k}{ Integer, number of suboptimal paths to identify. }
  \item{ncore}{ Number of CPU cores used to do the calculation. 
       By default (NULL), use all detected CPU cores. }
  \item{\dots}{ Additional arguments passed to igraph function 
       \code{\link[igraph:get.shortest.paths]{get.shortest.paths}} (in 
       the function \code{cnapath}), passed to \code{summary.cnapath} 
       (in \code{print.cnapath}), as additional paths for comparison 
       (in \code{summary.cnapath}), or passed to the function \code{colorRamp} 
       (in \code{view.cnapath}). }
  \item{object}{ A \sQuote{cnapath} class of object obtained from 
       \code{\link{cnapath}}. Multiple \sQuote{pa} input is allowed for 
       comparative statistical analysis in \code{\link{summary.cnapath}}. }
  \item{pdb}{ A \sQuote{pdb} class of object obtained from \code{\link{read.pdb}} 
       and is used as the reference for node residue ids or for molecular 
       visulaization in VMD. }
  \item{label}{ Character, labels for paths identified from different networks. }
  \item{col}{ colors for plotting statistical resutls for paths identified 
       from different networks. }
  \item{plot}{ logical, if TRUE path length distribution and node degeneracy will be plotted. }
  \item{concise}{ logical, if TRUE only \sQuote{on path} residues will be displayed. }
  \item{cutoff}{ numeric, degeneracy cutoff for displaying nodes on paths. }
  \item{normalize}{ logical, if TRUE node degeneracy is defined by the percentage 
    of total number of paths. }
  \item{x}{ A \sQuote{cnapath} class of object obtained from 
    \code{\link{cnapath}}. }
  \item{out.prefix}{ prefix for the names of output files, \sQuote{view.vmd}
       and \sQuote{view.pdb}. }
  \item{spline}{ logical, if TRUE all paths are displayed as spline curves. }
  \item{colors}{ character vector or integer scalar, define path colors. If a
     character vector, passed to \code{colorRamp} function to generate the 
     color scales. If an integer, color all paths the same way with VMD color 
     ID equal to the integer. }
  \item{launch}{ logical, if TRUE VMD will be launched. }
}
\value{
  The function \code{\link{cnapath}} returns a \sQuote{cnapath} class of list containing following three components:
  \item{path}{ a list object containing all identified suboptimal paths.  
      Each entry of the list is a sequence of node ids for the path. }
  \item{epath}{ a list object containing all identified suboptimal paths. 
      Each entry of the list is a sequence of edge ids for the path. }
  \item{dist}{ a numeric vector of all path lengths. }

  The function \code{\link{summary.cnapath}} with return a matrix of (normalized)
  node degeneracy for \sQuote{on path} residues. 
}
\references{
  Yen, J.Y. (1971) \emph{Management Science} \bold{17}, 712--716.
}
\author{ Xin-Qiu Yao }
\seealso{ \code{\link{cna}}, \code{\link{cna.dccm}}, 
   \code{\link[igraph:get.shortest.paths]{get.shortest.paths}}. }
\examples{
\donttest{

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
#view.cnapath(pa1, pdb, launch = TRUE)
#view.cnapath(pa1, pdb, colors = 7, launch = TRUE)
#view.cnapath(pa1, pdb, spline=TRUE, colors=c("pink", "red"), launch = TRUE)
#pdb2 <- read.pdb("1tag")
#pdb2 <- trim.pdb(pdb2, atom.select(pdb2, chain="A", resno=npdbs$resno[2, gaps.res$f.inds]))
#view.cnapath(pa2, pdb2, launch = TRUE)

detach(transducin)
}
}
\keyword{ utilities }
