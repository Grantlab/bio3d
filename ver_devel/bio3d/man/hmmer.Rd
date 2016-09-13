\name{hmmer}
\alias{hmmer}
\title{ HMMER Sequence Search }
\description{
  Perform a HMMER search against the PDB, NR, swissprot or other sequence
  and structure databases.
}
\usage{
hmmer(seq, type="phmmer", db = NULL, verbose = TRUE, timeout = 90)
}
\arguments{
  \item{seq}{ a multi-element character vector containing the query
    sequence. Alternatively a \sQuote{fasta} object as obtained
    from functions \code{get.seq} or \code{read.fasta} can be
    provided. }
  \item{type}{ character string specifying the \sQuote{HMMER} job
    type. Current options are \sQuote{phmmer}, \sQuote{hmmscan},
    \sQuote{hmmsearch}, and \sQuote{jackhmmer}. }
  \item{db}{ character string specifying the database to search.
    Current options are \sQuote{pdb}, \sQuote{nr}, \sQuote{swissprot},
    \sQuote{pfam}, etc. See \sQuote{details} for a complete list. }
  \item{verbose}{ logical, if TRUE details of the download process
    is printed. }
  \item{timeout}{ integer specifying the number of seconds to wait for
    the blast reply before a time out occurs. }
}
\details{
  This function employs direct HTTP-encoded requests to the HMMER web
  server. HMMER can be used to search sequence databases for homologous
  protein sequences. The HMMER server implements methods using
  probabilistic models called profile hidden Markov models (profile
  HMMs).

  There are currently four types of HMMER search to perform:
  
  - \sQuote{phmmer}: protein sequence vs protein sequence database.\cr
  (input argument \code{seq} must be a sequence). 
  
  Allowed options for \code{type} includes: 
  \sQuote{env_nr}, \sQuote{nr}, \sQuote{refseq}, \sQuote{pdb},
  \sQuote{rp15}, \sQuote{rp35}, \sQuote{rp55}, \sQuote{rp75},
  \sQuote{swissprot}, \sQuote{unimes}, \sQuote{uniprotkb},
  \sQuote{uniprotrefprot}, \sQuote{pfamseq}. 

  
  - \sQuote{hmmscan}: protein sequence vs profile-HMM database.\cr
  (input argument \code{seq} must be a sequence). 
  
  Allowed options for \code{type} includes:
  \sQuote{pfam}, \sQuote{gene3d},  \sQuote{superfamily}, \sQuote{tigrfam}.
  

  - \sQuote{hmmsearch}: protein alignment/profile-HMM vs protein sequence
  database.\cr
  (input argument \code{seq} must be an alignment). 
  
  Allowed options for \code{type} includes:
  \sQuote{pdb}, \sQuote{swissprot}.
  

  - \sQuote{jackhmmer}: iterative search vs protein sequence database.\cr
  (input argument \code{seq} must be an alignment).
  \sQuote{jackhmmer} functionality incomplete!!
  
  Allowed options for \code{type} includes:
  \sQuote{env_nr}, \sQuote{nr}, \sQuote{refseq}, \sQuote{pdb},
  \sQuote{rp15}, \sQuote{rp35}, \sQuote{rp55}, \sQuote{rp75},
  \sQuote{swissprot}, \sQuote{unimes}, \sQuote{uniprotkb},
  \sQuote{uniprotrefprot}, \sQuote{pfamseq}. 
  
  
  More information can be found at the HMMER website:\cr
  \url{http://hmmer.janelia.org}
}
\value{
  A data frame with multiple components depending on the selected job
  \sQuote{type}. Frequently reported fields include:

  \item{name}{ a character vector containing the name of the target. }
  \item{acc}{ a character vector containing the accession identifier of
    the target. }
  \item{acc2}{ a character vector containing secondary accession of the
    target. }
  \item{id}{ a character vector containing Identifier of the target }
  \item{desc}{ a character vector containing entry description. }
  \item{score}{ a numeric vector containing bit score of the sequence (all
    domains, without correction). }
  \item{pvalue}{ a numeric vector containing the P-value of the score. }
  \item{evalue}{ a numeric vector containing the E-value of the score. }
  \item{nregions}{ a numeric vector containing Number of regions
    evaluated. } 
  \item{nenvelopes}{ a numeric vector containing the number of envelopes
    handed over for domain definition, null2, alignment, and scoring. }
  \item{ndom}{ a numeric vector containing the total number of domains
    identified in this sequence. } 
  \item{nreported}{ a numeric vector containing the number of domains
    satisfying reporting thresholding. } 
  \item{nincluded}{ a numeric vector containing the number of domains
    satisfying inclusion thresholding. } 
  \item{taxid}{ a character vector containing The NCBI taxonomy identifier
    of the target (if applicable). } 
  \item{species}{ a character vector containing the species name. }
  \item{kg}{ a character vector containing the kingdom of life that the
    target belongs to - based on placing in the NCBI taxonomy tree. }

  More details can be found at the HMMER website:\cr
  \url{http://www.ebi.ac.uk/Tools/hmmer/help/api}
}
\note{
  Note that the chained \sQuote{pdbs} HMMER field (used for redundant
  PDBs) is included directly into the result list (applies only when
  \code{db='pdb'}). In this case, the \sQuote{name} component of the
  target contains the parent (non redundant) entry, and the \sQuote{acc}
  component the chained PDB identifiers. The search results will therefore
  provide duplicated PDB identifiers for component \code{$name}, while
  \code{$acc} should be unique. 
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.

  Finn, R.D. et al. (2011) \emph{Nucl. Acids Res.} \bold{39}, 29--37.
  Eddy, S.R. (2011) \emph{PLoS Comput Biol} \bold{7}(10): e1002195.
  
  See also the \sQuote{HMMER} website:\cr
  \url{http://hmmer.org}
}
\author{ Lars Skjaerven }
\note{
  Online access is required to query HMMER services.
}
\seealso{
  \code{\link{seqaln}}, \code{\link{get.seq}}
  \code{\link{blast.pdb}}, \code{\link{pfam}}
  \code{\link{uniprot}}
}
\examples{
\dontrun{
# HMMER server connection required - testing excluded

##- PHMMER
seq <- get.seq("2abl_A", outfile=tempfile())
res <- hmmer(seq, db="pdb")

##- HMMSCAN
fam <- hmmer(seq, type="hmmscan", db="pfam")
pfam.aln <- pfam(fam$acc[1])

##- HMMSEARCH
hmm <- hmmer(pfam.aln, type="hmmsearch", db="pdb")
unique(hmm$species)
hmm$acc

}
}
\keyword{ utilities }
