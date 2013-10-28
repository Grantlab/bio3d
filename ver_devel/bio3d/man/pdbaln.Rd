\name{pdbaln}
\alias{pdbaln}
\title{ Sequence Alignment of PDB Files }
\description{
  Create multiple sequences alignments from a list of PDB files
  returning aligned sequence and structure records.
}
\usage{
pdbaln(files, fit = TRUE, pqr = FALSE, ncore = 1, nseg.scale = 1,  ...) 
}
\arguments{
  \item{files}{ a character vector of PDB file names. }
  \item{fit}{ logical, if TRUE coordinate superposition is performed on
    the input structures. }
  \item{pqr}{ logical, if TRUE the input structures are assumed to be in
  PQR format. }
  \item{ncore }{ number of CPU cores used to do the calculation.
    \code{ncore>1} requires multicore package installed. }
  \item{nseg.scale }{ split input data into specified number of segments
    prior to running multiple core calculation. See \code{\link{fit.xyz}}. }
  \item{\dots}{ extra arguments passed to \code{seqaln} function. }
}
\details{
  This wrapper function calls the underlying functions \code{read.pdb},
  \code{pdbseq}, \code{seqaln} and \code{read.fasta.pdb} returning a
  list of class \code{"3dalign"} similar to that returned by
  \code{read.fasta.pdb}.

  As these steps are often error prone it is recomended for most cases that
  the individual underlying functions are called in sequence with checks
  made on the valadity of their respective outputs to ensure sensible
  results.
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
\note{ See recommendation in details section above. }
\seealso{
  \code{\link{read.pdb}}, \code{\link{pdbseq}}, \code{\link{seqaln}},
  \code{\link{read.fasta}},\code{\link{read.fasta.pdb}},
  \code{\link{core.find}}, \code{\link{fit.xyz}}, \code{\link{read.all}}
}
\examples{
\dontrun{
#files <- get.pdb(c("4q21","5p21"), URLonly=TRUE)
files <- get.pdb(c("4q21","5p21"), path=tempdir(), overwrite=TRUE)
pdbaln(files)
}
}
\keyword{ utilities }
