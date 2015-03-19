\name{blast.pdb}
\alias{blast.pdb}
\alias{get.blast}
\alias{plot.blast}
\title{ NCBI BLAST Sequence Search and Summary Plot of Hit Statistics}
\description{
  Run NCBI blastp, on a given sequence, against the PDB, NR and
  swissprot sequence databases. Produce plots that facilitate hit selection
  from the match statistics of a BLAST result.
}
\usage{
blast.pdb(seq, database = "pdb", time.out = NULL, chain.single=TRUE)

get.blast(urlget, time.out = NULL, chain.single=TRUE)

\method{plot}{blast}(x, cutoff = NULL, cut.seed=NULL, cluster=TRUE, mar=c(2, 5, 1, 1), cex=1.5, ...)
}
\arguments{
  \item{seq}{ a single element or multi-element character vector
    containing the query sequence. Alternatively a \sQuote{fasta}
    object from function \code{get.seq} or \sQuote{pdb} object from 
    function \code{read.pdb} can be provided. }
  \item{database}{ a single element character vector specifying the
    database against which to search. Current options are \sQuote{pdb},
    \sQuote{nr} and \sQuote{swissprot}. }
  \item{time.out}{ integer specifying the number of seconds to wait for
    the blast reply before a time out occurs. }
  \item{urlget}{ the URL to retrieve BLAST results; Usually it is returned
    by blast.pdb if time.out is set and met. }
  \item{chain.single}{ logical, if TRUE double NCBI character PDB database  
    chain identifiers are simplified to lowercase '1WF4_GG' > '1WF4_g'. If 
    FALSE no conversion to match RCSB PDB files is performed. } 
  \item{x}{ BLAST results as obtained from the function
    \code{\link{blast.pdb}}. }
  \item{cutoff}{ A numeric cutoff value, in terms of minus the log of
    the evalue, for returned hits.  If null then the function will
  try to find a suitable cutoff near \sQuote{cut.seed} which can be used
  as an initial guide (see below). }
  \item{cut.seed}{ A numeric seed cutoff value, used for initial
    cutoff estimation. If null then a seed position is set to the point of 
    largest drop-off in normalized scores (i.e. the biggest jump in 
    E-values). }
  \item{cluster}{ Logical, if TRUE (and \sQuote{cutoff} is null) a clustering 
    of normalized scores is performed to partition hits in groups by 
    similarity to query. If FALSE the partition point is set to the point of largest drop-off in normalized scores. }
  \item{mar}{ A numerical vector of the form c(bottom, left, top, right)
    which gives the number of lines of margin to be specified on the
    four sides of the plot.}
  \item{cex}{ a numerical single element vector giving the amount by
    which plot labels should be magnified relative to the default. }
  \item{\dots}{ extra plotting arguments. }    
}
\details{
  The \code{blast.pdb} function employs direct HTTP-encoded requests to the NCBI web
  server to run BLASTP, the protein search algorithm of the BLAST
  software package.

  BLAST, currently the most popular pairwise sequence comparison algorithm 
  for database searching, performs gapped local alignments via a heuristic 
  strategy: it identifies short nearly exact matches or hits, bidirectionally 
  extends non-overlapping hits resulting in ungapped extended hits or 
  high-scoring segment pairs(HSPs), and finally extends the highest scoring 
  HSP in both directions via a gapped alignment (Altschul et al., 1997)

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

  Examining plots of BLAST alignment lengths, scores, E-values and normalized
  scores (-log(E-Value) from the \code{blast.pdb} function can aid in the
  identification sensible hit similarity thresholds. This is facilitated by 
  the \code{plot.blast} function.

  If a \sQuote{cutoff} value is not supplied then a basic hierarchical
  clustering of normalized scores is performed with initial group partitioning
  implemented at a hopefully sensible point in the vicinity of
  \sQuote{h=cut.seed}. Inspection of the resultant plot can then be use
  to refine the value of \sQuote{cut.seed} or indeed
  \sQuote{cutoff}. As the \sQuote{cutoff} value can vary depending on
  the desired application and indeed the properties of the system under
  study it is envisaged that \sQuote{plot.blast} will be called multiple
  times to aid selection of a suitable \sQuote{cutoff} value. See the
  examples below for further details.
}
\value{
  The function \code{blast.pdb} returns a list with the first eight components 
  below.  The function \code{plot.blast} produces a plot on the active graphics 
  device and returns a three component list object with \code{hits}, \code{pdb.id} and \code{gi.id} see below:
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
  \item{url }{ a single element character vector with the NCBI result URL and 
    RID code. This can be passed to the get.blast function. }  
  \item{hits}{ an ordered matrix detailing the subset of hits with a normalized
    score above the chosen cutoff. Database identifiers are listed along
    with their cluster group number. }
  \item{pdb.id}{ a character vector containing the PDB database identifier of
    each hit above the chosen threshold. }
  \item{gi.id}{ a character vector containing the gi database identifier of
    each hit above the chosen threshold. }
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
\seealso{ \code{\link{plot.blast}}, \code{\link{hmmer}}, \code{\link{seqaln}}, \code{\link{get.pdb}} }
\examples{
\dontrun{
pdb <- read.pdb("4q21")
blast <- blast.pdb( pdbseq(pdb) )

head(blast$hit.tbl)
top.hits <- plot(blast)
head(top.hits$hits)

## Use 'get.blast()' to retrieve results at a later time.
#x <- get.blast(blast$url)
#head(x$hit.tbl)

# Examine and download 'best' hits
top.hits <- plot.blast(blast, cutoff=188)
head(top.hits$hits)
#get.pdb(top.hits)
}
}
\keyword{ utilities }
\keyword{ hplot }

