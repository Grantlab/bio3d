\name{get.seq}
\alias{get.seq}
\title{ Download FASTA Sequence Files }
\description{
  Downloads FASTA sequence files from the NR, or SWISSPROT/UNIPROT
  databases.
}
\usage{
get.seq(ids, outfile = "seqs.fasta", db = "nr", verbose = FALSE)
}
\arguments{
  \item{ids}{ A character vector of one or more appropriate database
    codes/identifiers of the files to be downloaded. }
  \item{outfile}{ A single element character vector specifying the name
    of the local file to which sequences will be written. }
  \item{db}{ A single element character vector specifying the database
    from which sequences are to be obtained. }
  \item{verbose}{ logical, if TRUE URL details of the download process
    are printed, else a progress bar is displayed. }
}
\details{
  This is a basic function to automate sequence file download from the
  NR and SWISSPROT/UNIPROT databases.
}
\value{
  If all files are successfully downloaded a list object with two
  components is returned:
  \item{ali }{ an alignment character matrix with a row per sequence and
    a column per equivalent aminoacid/nucleotide. }
  \item{ids }{ sequence names as identifiers.}
  This is similar to that returned by \code{\link{read.fasta}}. However,
  if some files were not successfully downloaded then a vector detailing
  which ids were not found is returned.
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Barry Grant }
\note{  
  For a description of FASTA format see:
  \url{http://www.ncbi.nlm.nih.gov/BLAST/blastcgihelp.shtml}.
  When reading alignment files, the dash \sQuote{-} is interpreted as
  the gap character.
}
\seealso{ \code{\link{blast.pdb}}, \code{\link{read.fasta}},
  \code{\link{read.fasta.pdb}}, \code{\link{get.pdb}}
}
\examples{
\dontrun{
## Sequence identifiers (GI or PDB codes e.g. from blast.pdb etc.)
get.seq( c("P01112", "Q61411", "P20171") )

#aa <-get.seq( c("4q21", "5p21") )
#aa$id
#aa$ali
}
}
\keyword{ IO }
\keyword{ utilities }
