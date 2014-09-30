\name{straln}
\alias{straln}
\title{
  Multiple Structure Alignment
}
\description{
  This function applies \code{straln.pair} and does tree-guided 
  multiple structure alignment.
}
\usage{
  straln(pdbs, prefix = "", outfile = "straln.fa", ret.hc = FALSE, 
         ncore = NULL, ...)
}
\arguments{
  \item{pdbs}{ A \sQuote{3dalign} object as obtained from calling \code{read.fasta.pdb}. }
  \item{prefix}{ Prefix of the file paths contained in \sQuote{pdbs}. }
  \item{outfile}{ Output alignment file name. }
  \item{ret.hc}{ logical, if TRUE returns both alignment and information about clustering. }
  \item{ncore}{ number of CPU cores used to do the calculation. By default, all cores will be used}
  \item{\dots}{ extra arguments passed to \code{straln.pair} function. }
}
\value{
  if \sQuote{ret.hc} is FALSE, returns an alignment object similar to \code{seqaln}.
  if \sQuote{ret.hc} is TRUE, returns:
  \item{aln}{ An alignment object. }
  \item{hc}{ Tree built by calling \code{hclust}. }
  \item{distmat}{ The distance matrix used for clustering. }
}
\note{ 
  A couple things could be improved such as the support of different tree building algorithms including
  hierarchical clustering, Neighbor-Joining, UPGMA, etc.
}  
\seealso{
   \code{\link{straln.pair}}, \code{\link{seqaln}}, \code{\link{seqaln.pair}}, \code{\link{struct.aln}}
}
\examples{

\dontrun{

# Pair-wise structure alignment
#files <- get.pdb(c("1as0","5p21"), path=tempdir(), overwrite=TRUE)
#pdbs <- pdbaln(files)
#pdbs$xyz <- pdbfit(pdbs)
#pdb <- pdbs2pdb(pdbs, inds=c(1,2))
#aln <- straln.pair(pdb[[1]], pdb[[2]])

# Profile alignment

}
}
\keyword{analysis}
