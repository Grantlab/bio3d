\name{nma}
\alias{nma}
\title{ Normal Mode Analysis }
\description{
  Perform normal mode analysis (NMA) on either a single or an ensemble
  of protein structures.
}
\usage{
nma(...)
}
\arguments{
  \item{\dots}{ arguments passed to the methods
    \code{\link{nma.pdb}}, or \code{\link{nma.pdbs}}.
    
    For function \code{\link{nma.pdb}} this will include an object of
    class \code{pdb} as obtained from function \code{\link{read.pdb}}.
    
    For function \code{\link{nma.pdbs}} an object of class \code{pdbs} as
    obtained from function \code{\link{pdbaln}} or
    \code{\link{read.fasta.pdb}}. } 
}
\details{
  Normal mode analysis (NMA) is a computational approach for studying and
  characterizing protein flexibility. Current functionality entails normal
  modes calculation on either a single protein structure or an ensemble
  of aligned protein structures. 

  This generic \code{\link{nma}} function calls the corresponding
  methods for the actual calculation, which is determined by the class
  of the input argument:

  Function \code{\link{nma.pdb}} will be used when the input argument is
  of class \code{pdb}. The function calculates the normal modes of a
  C-alpha model of a protein structure.

  Function \code{\link{nma.pdbs}} will be used when the input argument is
  of class \code{pdbs}. The function will perform normal mode analysis
  of each PDB structure stored in the \code{pdbs} object
  (\sQuote{ensemble NMA}).
  
  See documentation and examples for each corresponding function for
  more details. 
}
\references{
  Skjaerven, L. et al. (2014) \emph{BMC Bioinformatics} \bold{15}, 399.
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Lars Skjaerven }
\seealso{
  \code{\link{nma.pdb}}, \code{\link{nma.pdbs}},
  \code{\link{pca}}.
}
\examples{

##- Singe structure NMA
## Fetch stucture
pdb <- read.pdb( system.file("examples/1hel.pdb", package="bio3d") )

## Calculate normal modes
modes <- nma(pdb)

## Print modes
print(modes)

## Plot modes
plot(modes)

## Visualize modes
#m7 <- mktrj.nma(modes, mode=7, file="mode_7.pdb")

\donttest{
## Needs MUSCLE installed - testing excluded

##- Ensemble NMA
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

}
}
}
\keyword{ analysis }
