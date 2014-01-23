\name{nma.pdbs}
\alias{nma.pdbs}
\alias{print.enma}
\title{ Ensemble Normal Mode Analysis }
\description{
  Perform NMA on an ensemble of aligned protein structures.
}
\usage{
nma.pdbs(pdbs, fit=TRUE, full=FALSE, rm.gaps=TRUE, defa=FALSE, 
         outpath="pdbs_nma", \dots)

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
  \item{rm.gaps}{ logical, if TRUE obtain the hessian matrices for only
    atoms in the aligned positions (non-gap positions in all aligned
    structures). Thus, gap positions are removed from output. }
  \item{defa}{ logical, if TRUE perform deformation analysis on the
    first 5 modes. }
  \item{outpath}{ character string specifing the output directory to
    which the PDB structures should be written. }
  \item{x}{ an \code{enma} object obtained from \code{\link{nma.pdbs}}. }
  \item{...}{ additional arguments to \code{nma}, \code{aa2mass}, and
    \code{print}. }
}
\details{
  This function performs normal mode analysis (NMA) on a set of aligned
  protein structures obtained with function \code{read.fasta.pdb} or
  \code{pdbaln}. The main purpose is to provide aligned atomic
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
  \item{fluctuations }{ a numeric matrix with aligned atomic
    fluctuations. }
  \item{rmsip}{ a numeric matrix of RMSIP values between all pairs
    of mode subsets. }
  \item{deformations }{ a numeric matrix with aligned atomic
    deformation energies. }
  \item{U.subspace }{ a three-dimensional array with aligned
    eigenvectors  (corresponding to the subspace defined by the first 20
    non-trivial eigenvectors (\sQuote{U}) of the \sQuote{nma} object). }
  \item{full.nma }{ a list with a \code{nma} object for each input
    structure. }
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Lars Skjaerven }
\seealso{ \code{\link{nma}}, \code{\link{pdbaln}},
  \code{\link{read.fasta.pdb}}, \code{\link{rmsip}}  }
\examples{
## Fetch PDB files and split to chain A only PDB files
ids <- c("1a70_A", "1czp_A", "1frd_A", "1fxi_A", "1iue_A", "1pfd_A")
raw.files <- get.pdb(ids, path = "raw_pdbs")
files <- pdbsplit(raw.files, ids, path = "raw_pdbs/split_chain")

## Sequence Alignement
pdbs <- pdbaln(files)

## Normal mode analysis on aligned data
all.modes <- nma.pdbs(pdbs, rm.gaps=FALSE)

## Plot fluctuation data
plot.bio3d(all.modes$fluctuations[1,], type="o", xlab="Residue Position")
lines(all.modes$fluctuations[2,], type='o', col=2)
lines(all.modes$fluctuations[3,], type='o', col=3)

## Remove gaps from output
all.modes <- nma.pdbs(pdbs, rm.gaps=TRUE)

## Compare modes between two structures 
rmsip(all.modes$U.subs[,,1], all.modes$U.subs[,,2])

}
\keyword{ utilities }
