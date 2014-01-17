\name{dccm.mean}
\alias{dccm.mean}
\title{ Atomic Correlation Matrix Consensus Filtering }
\description{
Filter a three-dimensional correlation matrix, or DCCM, composed of 
results from multiple dccm calculations and return correlations 
present in at least a chosen subset of the calculated matrices.
}
\usage{
\method{dccm}{mean}(x, cutoff.sims = dim(x)[3], cutoff.cij = 0.4, ...)
}
\arguments{
  \item{x}{ A numeric array with 3 dimensions \sQuote{(NxNxZ)} containing  
  atomic correlation values for N residues or atoms from Z 
  calculations. } 
  \item{cutoff.sims}{ A single element numeric vector corresponding 
  to the minimum number of Z dimensions (i.e. dccm matrices) a 
  correlation between two residues must be present and whose absolute 
  value must be higher than \code{cutoff.cij}. If not, the matrix element 
  will be set to 0 and not used for the average calculation. }	
  \item{cutoff.cij}{ A single element numeric vector corresponding to 
  the minimum absolute value a correlation element must be equal to or 
  higher to be used for the average calculation. If lower, the element 
  will be set to 0 and not used for the average calculation. }
  \item{\dots}{ extra arguments that are currently ignored. }
}
\value{
Returns a numeric \sQuote{NxN} atomic correlation matrix of class \sQuote{dccm}.
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Guido Scarabelli & Barry Grant }
\seealso{
	\code{\link{dccm}}, \code{\link{dccm.nma}}, \code{\link{dccm.xyz}}, 
	\code{\link{plot.dccm}}
}
\examples{

\dontrun{
# Select Protein Kinase PDB IDs
ids <- c("4b7t_A", "2exm_A", "1opj_A", "4jaj_A", "1a9u_A", 
		 "1tki_A", "1csn_A", "1lp4_A")

# Download and split by chain ID
raw.files <- get.pdb(ids, path = "raw_pdbs") 
files <- pdbsplit(raw.files, ids)

# Alignment of structures
pdbs <- pdbaln(files) # Sequence identity
summary(c(seqidentity(pdbs)))

# NMA on all structures
modes <- nma.pdbs(pdbs, full = TRUE)

# Calculate correlation matrices for each structure
cij <- dccm(modes)

# Set DCCM plot panel names for combined figure
dimnames(cij$all.dccm) = list(NULL, NULL, ids) 
plot.dccm(cij$all.dccm)

# Filter to display only correlations present in all structures
cij.all <- dccm.mean(cij$all.dccm, cutoff.sims = 8, cutoff.cij = 0) 
plot.dccm(cij.all, main = "Consensus Residue Cross Correlation")
}
}
\keyword{analysis}

