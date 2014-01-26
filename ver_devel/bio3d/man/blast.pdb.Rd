\name{blast.pdb}
\alias{blast.pdb}
\title{ NCBI BLAST Sequence Search }
\description{
  Run NCBI blastp, on a given sequence, against the PDB, NR and
  swissprot sequence databases.
}
\usage{
blast.pdb(seq, database = "pdb")
}
\arguments{
  \item{seq}{ a single element or multi-element character vector
    containing the query sequence. Alternatively a \sQuote{fasta}
    object from function \code{get.seq} can be provided. }
  \item{database}{ a single element character vector specifying the
    database against which to search. Current options are \sQuote{pdb},
    \sQuote{nr} and \sQuote{swissprot}. }
}
\details{
  This function employs direct HTTP-encoded requests to the NCBI web
  server to run BLASTP, the protein search algorithm of the BLAST
  software package.

  BLAST, currently the fastest and most popular pairwise sequence
  comparison algorithm, performs gapped local alignments, through the
  implementation of a heuristic strategy: it identifies short nearly
  exact matches or hits, bidirectionally extends non-overlapping hits
  resulting in ungapped extended hits or high-scoring segment pairs
  (HSPs), and finally extends the highest scoring HSP in both directions
  via a gapped alignment (Altschul et al., 1997)

  For each pairwise alignment BLAST reports the raw score, bitscore
  and an E-value that assess the statistical significance of the
  raw score. Note that unlike the raw score E-values are normalized with
  respect to both the substitution matrix and the query and database lengths.

  Here we also return a corrected normalized score (mlog.evalue) that in
  our experience is easier to handle and store than conventional
  E-values. In practice, this score is equivalent to minus the natural
  log of the E-value. Note that, unlike the raw score, this score is
  independent of the substitution matrix and and the query and database
  lengths, and thus is comparable between BLASTP searches. 

}
\value{
  A list with seven components:
  \item{bitscore }{ a numeric vector containing the raw score for each
    alignment. }
  \item{evalue }{ a numeric vector containing the E-value of the raw
    score for each alignment. }
  \item{mlog.evalue }{ a numeric vector containing minus the natural log
    of the E-value. }
  \item{gi.id }{ a character vector containing the gi database identifier of
    each hit. }  
  \item{pdb.id }{ a character vector containing the PDB database identifier of
    each hit. }  
  \item{hit.tbl }{ a character matrix summarizing BLAST results for each
    reported hit, see below. }
  \item{raw }{ a data frame summarizing BLAST results, note multiple
    hits may appear in the same row. }  
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.

  \sQuote{BLAST} is the work of Altschul et al.:
  Altschul, S.F. et al. (1990) \emph{J. Mol. Biol.} \bold{215}, 403--410.
  
  Full details of the \sQuote{BLAST} algorithm, along with download and
  installation instructions can be obtained from:\cr
  \url{http://www.ncbi.nlm.nih.gov/BLAST/}.
}
\author{ Barry Grant }
\note{
  Online access is required to query NCBI blast services.
}
\seealso{ \code{\link{seqaln}} }
\examples{
pdb <- read.pdb("1bg2")
blast <- blast.pdb( pdbseq(pdb) )

head(blast$hit.tbl)
top.hits <- plot(blast)
head(top.hits$hits)
}
\keyword{ utilities }
