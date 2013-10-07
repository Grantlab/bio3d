\name{nma.pdbs}
\alias{nma.pdbs}
\title{ Ensemble Normal Mode Analysis }
\description{
  Perform NMA on an ensemble of aligned protein structures.
}
\usage{
nma.pdbs(aln, fit=FALSE, full=FALSE, trim.inds=NULL,
         outpath="pdbs_nma/", \dots)
}
\arguments{
  \item{aln}{ a numeric matrix of aligned C-alpha xyz Cartesian
    coordinates. For example an alignment data structure obtained with
    \code{\link{read.fasta.pdb}} or \code{\link{pdbaln}}. }
  \item{fit}{ logical, if TRUE coordinate superposition is performed
    prior to normal mode calculations. }
  \item{full}{ logical, if TRUE return the complete, full structure, 
    \sQuote{nma} objects. By default, only fluctuation data for aligned 
    positions in all structures is returned. }
  \item{trim.inds}{ a vector of alignment column indices that selects the 
    positions and corresponding coordinates in all aligned structures upon 
    which NMA will be performed. If provided, the PDB
    structures will be trimmed to obtain equal lengths. }
  \item{outpath}{ character string specifing the output directory to
    which the PDB structures should be written. }
  \item{...}{ additional arguments to \code{nma}. }
}
\details{
  This function performs normal mode analysis (NMA) on a set of aligned
  protein structures obtained with function
  \code{read.fasta.pdb} or \code{pdbaln}. The main purpose is to provide 
  aligned atomic fluctuations and mode vectors in an automated fashion. 

  NOTE: The normal modes are calculated on the full structures as provided 
  by object \sQuote{aln}. With the input argument \sQuote{full=TRUE} the 
  output \sQuote{modes.array} provides mode vectors in which unaligned atoms 
  are ommited. Thus, these mode vectors might not be orthogonal anymore.
  
  See examples for more details. 
 }
\value{
  Returns a list with the following components:
  \item{fluctuations }{ a numeric matrix with aligned atomic
    fluctuations. }
  \item{modes.array }{ a three-dimensional array with aligned mode
    vectors. }
  \item{all.modes }{ a list with a \code{nma} object for each input
    structure. }
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Lars Skjaerven }
\seealso{ \code{\link{nma}}, \code{\link{read.fasta.pdb}} }
\examples{
## Fetch PDB files and split to chain A only PDB files
ids <- c("1a70_A", "1czp_A", "1frd_A")
raw.files <- get.pdb(ids, path = "raw_pdbs")
files <- pdbsplit(raw.files, ids, path = "raw_pdbs/split_chain/")

## Sequence Alignement
pdbs <- pdbaln(files)

## Normal mode analysis on aligned data
all.modes <- nma.pdbs(pdbs, fit=TRUE)

## Plot fluctuation data
plot.bio3d(all.modes$fluctuations[1,], typ="o", xlab="Residue Position")
lines(all.modes$fluctuations[2,], type='o', col=2)
lines(all.modes$fluctuations[3,], type='o', col=3)

## Provide full mode objects
all.modes <- nma.pdbs(pdbs, fit=TRUE, full=TRUE)

## Compare modes between two structures 
rmsip(all.modes$modes.array[,,1], all.modes$modes.array[,,2])

}
\keyword{ utilities }
