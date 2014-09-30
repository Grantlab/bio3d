\name{hmmer}
\alias{hmmer}
\title{ HMMER Sequence Search }
\description{
  Run HMMER, against the PDB, NR, swissprot (and other) sequence
  databases. 
}
\usage{
hmmer(seq, type="phmmer", db = NULL, verbose = TRUE, timeout = 90)
}
\arguments{
  \item{seq}{ a single element or multi-element character vector
    containing the query sequence. Alternatively a \sQuote{fasta}
    object from function \code{get.seq} can be provided. }
  \item{type}{ for now: \sQuote{phmmer}, or \sQuote{hmmscan}.  
    Currently not available: \sQuote{hmmsearch}, \sQuote{jackhmmer} }
  \item{db}{ a single element character vector specifying the
    database against which to search. Current options are \sQuote{pdb},
    \sQuote{nr} and \sQuote{swissprot}. More here... }
  \item{verbose}{ logical, if TRUE details of the download process
    is printed. }
  \item{timeout}{ integer specifying the number of seconds to wait for
    the blast reply before a time out occurs. }
}
\details{
  This function employs direct HTTP-encoded requests to the HMMER web
  server ...

  Coppy from the HMMER website:
  HMMER is used for searching sequence databases for homologs of protein
  sequences, and for making protein sequence alignments. It implements
  methods using probabilistic models called profile hidden Markov
  models (profile HMMs).
}
\value{
  A data frame with the following components:
  this depends on the \sQuote{type}.
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.

  Finn, R.D. et al. (2011) \emph{Nucl. Acids Res.} \bold{39}, 29--37.
  Eddy, S.R. (2011) \emph{PLoS Comput Biol} \bold{7}(10): e1002195.
  
  See also the \sQuote{HMMER} website:\cr
  \url{http://hmmer.janelia.org}
}
\author{ Lars Skjaerven }
\note{
  Online access is required to query HMMER services.
}
\seealso{ \code{\link{seqaln}} }
\examples{
##- PHMMER
seq <- get.seq("2abl_A")
res <- hmmer(seq, db="pdb")

##- HMMSCAN
fam <- hmmer(seq, type="hmmscan", db="pfam")
pfam.aln <- pfam(fam$acc[1])

##- HMMSEARCH
hmm <- hmmer(pfam.aln, type="hmmsearch", db="pdb")
unique(hmm$species)
hmm$acc

}
\keyword{ utilities }
