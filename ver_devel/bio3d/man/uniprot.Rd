\name{uniprot}
\alias{uniprot}
\title{ Fetch UniProt Entry Data. }
\description{
  Fetch protein sequence and functional information from the UniProt
  database. 
}
\usage{
uniprot(accid)
}
\arguments{
  \item{accid}{ UniProt accession id. }
}
\details{
  This is a basic utility function for downloading information from the
  UniProt database. UniProt contains protein sequence and functional
  information.
}
\value{
  Returns a list object with the following components:
  \item{accession}{ a character vector with UniProt accession id's. }
  \item{name}{ abbreviated name. }
  \item{fullName}{ full recommended protein name. }
  \item{shortName }{ short protein name. }
  \item{sequence}{ protein sequence. }
  \item{gene}{ gene names. }
  \item{organism}{ organism. }
  \item{taxon}{ taxonomic lineage. }
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
  
  See also the UniProt web-site for more information:\cr
  \url{http://www.uniprot.org/}.
}
\author{ Lars Skjaerven }
\seealso{ \code{\link{blast.pdb}}, \code{\link{get.seq}} }
\examples{
\dontrun{
# UNIPROT server connection required - testing excluded

prot <- uniprot('PH4H_HUMAN')
prot$fullName
prot$sequence
}
}
\keyword{ utilities }
