\name{read.fasta.pdb}
\alias{read.fasta.pdb}
\title{ Read Aligned Structure Data }
\description{
  Read aligned PDB structures and store their C-alpha atom data, including xyz
  coordinates, residue numbers, residue type and B-factors.
}
\usage{
read.fasta.pdb(aln, pdb.path = "", pdbext = "", 
             ncore = 1, nseg.scale = 1, ...)
}
\arguments{
  \item{aln}{ an alignment data structure obtained with
    \code{\link{read.fasta}}. }
  \item{pdb.path}{ path to PDB files. }
  \item{pdbext}{ the file name extention of the PDB files. }
  \item{ncore }{ number of CPU cores used to do the calculation.
    \code{ncore>1} requires multicore package installed. }
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
  Returns a list of class \code{"3dalign"} with the following five
  components:
  \item{xyz}{numeric matrix of aligned C-alpha coordinates.}
  \item{resno}{character matrix of aligned residue numbers.}
  \item{b}{numeric matrix of aligned B-factor values.}
  \item{chain}{character matrix of aligned chain identifiers.}
  \item{id}{character vector of PDB sequence/structure names.}
  \item{ali}{character matrix of aligned sequences.}
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
# Read sequence alignment
file <- system.file("examples/kif1a.fa",package="bio3d")
aln  <- read.fasta(file)

# Read aligned PDBs
pdb.path = system.file("examples",package="bio3d")
pdbs <- read.fasta.pdb(aln, pdb.path = pdb.path, pdbext = ".ent")

# Structure/sequence names/ids
pdbs$id

# Alignment positions 335 to 339
pdbs$ali[,335:339]
pdbs$resno[,335:339]
pdbs$b[,335:339]

# Alignment C-alpha coordinates for first structure
pdbs$xyz[1,]

}
\keyword{ IO }
