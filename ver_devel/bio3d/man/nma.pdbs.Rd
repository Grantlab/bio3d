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
    \code{\link{read.fasta.pdb}}. }
  \item{fit}{ logical, if TRUE coordinate superposition is performed
    prior to normal mode calculations. }
  \item{full}{ logical, if TRUE return the complete \sQuote{nma}
    objects. By default, only aligned fluctuation data is provided. }
  \item{trim.inds}{ a vector of indices that selects the elements upon
    which the calculation should be based. If provided, the PDB
    structures will be trimmed to obtain equal lengths. }
  \item{outpath}{ character string specifing the output directory to
    which the PDB structures should be written. }
  \item{...}{ additional arguments to \code{nma}. }
}
\details{
  This function performs normal mode analysis (NMA) on a set of aligned
  protein structures obtained with function
  \code{read.fasta.pdb}. The main purpose is to provide aligned
  atomic fluctuations and mode vectors in an automated fashion. 

  The normal modes are calculated on the full structures as provided by 
  object \sQuote{aln}. When \sQuote{full=TRUE} \sQuote{modes.array}
  provides mode vectors in which unaligned atoms are ommited. Thus,
  these mode vectors might not be orthogonal anymore.
  
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
## pdbIDs-chainIDs
ids <- c("1a70_a", "1czp_a", "1frd_a")

pdb.ids <- substr(ids, 1, 4)
cha.ids <- toupper(substr(ids, 6, 6))
ids     <- paste(pdb.ids, cha.ids, sep="_")

## Fetch and split PDB files
raw.files <- get.pdb(pdb.ids, path = "raw_pdbs")
files <- pdbsplit(raw.files, ids, path = "raw_pdbs/split_chain/")

## Sequence alignement
aln <- pdbaln(files)
pdbs <- read.fasta.pdb(aln, "", "", het2atom = T)

## Normal mode analysis on aligned data
all.modes <- nma.pdbs(pdbs, fit=TRUE)

## Plot fluctuation data
plot.bio3d(all.modes$fluctuations[1,])
lines(all.modes$fluctuations[2,], type='b', col=2)
lines(all.modes$fluctuations[3,], type='b', col=3)

## Provide full mode objects
all.modes <- nma.pdbs(pdbs, fit=TRUE, full=TRUE)

## Compare modes between two structures 
rmsip(all.modes$modes.array[,,1], all.modes$modes.array[,,2])

}
\keyword{ utilities }
