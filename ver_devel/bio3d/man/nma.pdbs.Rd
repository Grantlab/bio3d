\name{nma.pdbs}
\alias{nma.pdbs}
\alias{print.enma}
\title{ Ensemble Normal Mode Analysis }
\description{
  Perform normal mode analysis (NMA) on an ensemble of aligned protein
  structures.
}
\usage{
\method{nma}{pdbs}(pdbs, fit = TRUE, full = FALSE, subspace = NULL,
         rm.gaps = TRUE, varweight=FALSE, 
         outpath = NULL, ncore = 1, progress = NULL, \dots)

\method{print}{enma}(x, \dots)
}
\arguments{
  \item{pdbs}{ a numeric matrix of aligned C-alpha xyz Cartesian
    coordinates. For example an alignment data structure obtained with
    \code{\link{read.fasta.pdb}} or \code{\link{pdbaln}}. }
  \item{fit}{ logical, if TRUE coordinate superposition is performed
    prior to normal mode calculations. }
  \item{full}{ logical, if TRUE return the complete, full structure, 
    \sQuote{nma} objects. }
  \item{subspace}{ number of eigenvectors to store for further
    analysis. }
  \item{rm.gaps}{ logical, if TRUE obtain the hessian matrices for only
    atoms in the aligned positions (non-gap positions in all aligned
    structures). Thus, gap positions are removed from output. }
  \item{varweight}{ logical, if TRUE perform weighing of the pair force
    constants. Alternatively, provide a NxN matrix containing the
    weights. See function \code{\link{var.xyz}}. }
  \item{outpath}{ character string specifing the output directory to
    which the PDB structures should be written. }
  \item{ncore }{ number of CPU cores used to do the calculation.
    \code{ncore>1} requires package \sQuote{parallel} installed. }
  
  \item{x}{ an \code{enma} object obtained from \code{\link{nma.pdbs}}. }
  \item{progress}{ progress bar for use with shiny web app. }
  \item{...}{ additional arguments to \code{nma}, \code{\link{aa2mass}}, and
    \code{\link{print}}. }
}
\details{
  This function performs normal mode analysis (NMA) on a set of aligned
  protein structures obtained with function \code{\link{read.fasta.pdb}} or
  \code{\link{pdbaln}}. The main purpose is to provide aligned atomic
  fluctuations and mode vectors in an automated fashion. 

  The normal modes are calculated on the full structures as provided 
  by object \sQuote{pdbs}. With the input argument \sQuote{full=TRUE}
  the full \sQuote{nma} objects are returned together with output
  \sQuote{U.subs} providing the aligned mode vectors. When
  \sQuote{rm.gaps=TRUE} the unaligned atoms are ommited from
  output. With default arguments \sQuote{rmsip} provides RMSIP
  values for all pairwise structures.
  
  See examples for more details. 
 }
\value{
  Returns an \sQuote{enma} object with the following components:
  \item{fluctuations }{ a numeric matrix containing aligned atomic
    fluctuations with one row per input structure. }
  \item{rmsip}{ a numeric matrix of pair wise RMSIP values (only the ten
    lowest frequency modes are included in the calculation). }
  \item{U.subspace }{ a three-dimensional array with aligned
    eigenvectors  (corresponding to the subspace defined by the first N
    non-trivial eigenvectors (\sQuote{U}) of the \sQuote{nma} object). }
  \item{L}{ numeric matrix containing the raw eigenvalues with one row
    per input structure. }
  \item{xyz }{ an object of class \sQuote{xyz} containing the Cartesian
    coordinates in which the calculation was performed. Coordinates are
    superimposed to the first structure of the \code{pdbs} object when
    \sQuote{fit=TRUE}. }   
  \item{full.nma }{ a list with a \code{nma} object for each input
    structure. }
}
\references{
  Skjaerven, L. et al. (2014) \emph{BMC Bioinformatics} \bold{15}, 399.
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Lars Skjaerven }
\seealso{
  For normal mode analysis on single structure PDB:
  \code{\link{nma.pdb}}
  
  For the analysis of the resulting \sQuote{eNMA} object:
  \code{\link{mktrj.enma}}, \code{\link{dccm.enma}},
  \code{\link{plot.enma}}, \code{\link{cov.enma}}.
    
  Similarity measures:
  \code{\link{sip}}, \code{\link{covsoverlap}},
  \code{\link{bhattacharyya}}, \code{\link{rmsip}}.
  
  Related functionality:
  \code{\link{pdbaln}}, \code{\link{read.fasta.pdb}}.
}
\examples{
\donttest{
# Needs MUSCLE installed - testing excluded

if(check.utility("muscle")) {

## Fetch PDB files and split to chain A only PDB files
ids <- c("1a70_A", "1czp_A", "1frd_A", "1fxi_A", "1iue_A", "1pfd_A")
files <- get.pdb(ids, split = TRUE, path = tempdir())

## Sequence Alignement
pdbs <- pdbaln(files, outfile = tempfile())

## Normal mode analysis on aligned data
modes <- nma(pdbs, rm.gaps=FALSE)

## Plot fluctuation data
plot(modes, pdbs=pdbs)

## Cluster on Fluctuation similariy
sip <- sip(modes)
hc <- hclust(dist(sip))
col <- cutree(hc, k=3)

## Plot fluctuation data
plot(modes, pdbs=pdbs, col=col)

## Remove gaps from output
modes <- nma(pdbs, rm.gaps=TRUE)

## RMSIP is pre-calculated
heatmap(1-modes$rmsip)

## Bhattacharyya coefficient
bc <- bhattacharyya(modes)
heatmap(1-bc)

}
}
}
\keyword{ analysis }
