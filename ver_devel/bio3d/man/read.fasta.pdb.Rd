\name{read.fasta.pdb}
\alias{read.fasta.pdb}
\title{ Read Aligned Structure Data }
\description{
  Read aligned PDB structures and store their C-alpha atom data, including xyz
  coordinates, residue numbers, residue type and B-factors.
}
\usage{
read.fasta.pdb(aln, prefix = "", pdbext = "", fix.ali = FALSE,
             ncore = 1, nseg.scale = 1, ...)
}
\arguments{
  \item{aln}{ an alignment data structure obtained with
    \code{\link{read.fasta}}. }
  \item{prefix}{ prefix to aln$id to locate PDB files. }
  \item{pdbext}{ the file name extention of the PDB files. }
  \item{fix.ali}{ logical, if TRUE check consistence between \code{$ali} 
    and \code{$resno}, and correct \code{$ali} if they don't match. }
  \item{pdblist}{ an optional list of \code{pdb} objects with sequence
    corresponding to the alignments in \code{aln}. Primarily used
    through function \code{pdbaln} when the PDB objects already exists
    (avoids reading PDBs from file). }
  \item{ncore }{ number of CPU cores used to do the calculation.
    \code{ncore>1} requires package \sQuote{parallel} installed. }
  \item{nseg.scale }{ split input data into specified number of segments
    prior to running multiple core calculation. See \code{\link{fit.xyz}}. }
  \item{\dots}{ other parameters for \code{\link{read.pdb}}. }  
}
\details{
  The input \code{aln}, produced with \code{\link{read.fasta}}, must
  have identifers (i.e. sequence names) that match the PDB file
  names. For example the sequence corresponding to the structure
  \dQuote{1bg2.pdb} should have the identifer \sQuote{1bg2}. See
  examples below.
  
  Sequence miss-matches will generate errors.  Thus, care should be taken
  to ensure that the sequences in the alignment match the sequences in
  their associated PDB files.
}
\value{
  Returns a list of class \code{"pdbs"} with the following five
  components:
  \item{xyz}{numeric matrix of aligned C-alpha coordinates.}
  \item{resno}{character matrix of aligned residue numbers.}
  \item{b}{numeric matrix of aligned B-factor values.}
  \item{chain}{character matrix of aligned chain identifiers.}
  \item{id}{character vector of PDB sequence/structure names.}
  \item{ali}{character matrix of aligned sequences.}
  \item{resid}{character matrix of aligned 3-letter residue names.}
  \item{sse}{character matrix of aligned helix and strand secondary
    structure elements as defined in each PDB file.}  
  \item{call}{ the matched call. }
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Barry Grant }
\note{
   The sequence character \sQuote{X} is useful for masking unusual
   or unknown residues, as it can match any other residue type.
}
\seealso{
  \code{\link{read.fasta}}, \code{\link{read.pdb}},
  \code{\link{core.find}}, \code{\link{fit.xyz}}, \code{\link{read.all}} }
\examples{
\donttest{
# Redundant testing excluded

# Read sequence alignment
file <- system.file("examples/kif1a.fa",package="bio3d")
aln  <- read.fasta(file)

# Read aligned PDBs
pdbs <- read.fasta.pdb(aln)

# Structure/sequence names/ids
basename( pdbs$id )

# Alignment positions 335 to 339
pdbs$ali[,335:339]
pdbs$resid[,335:339]
pdbs$resno[,335:339]
pdbs$b[,335:339]

# Alignment C-alpha coordinates for these positions
pdbs$xyz[, atom2xyz(335:339)]

# See 'fit.xyz()' function for actual coordinate superposition
#  e.g. fit to first structure
# xyz <- fit.xyz(pdbs$xyz[1,], pdbs)
# xyz[, atom2xyz(335:339)]

}
}
\keyword{ IO }
