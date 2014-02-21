\name{straln.pair}
\alias{straln.pair}
\title{
  Pair-wise Structure Alignment with Dynamic Programming
}
\description{
  This function takes two \sQuote{3dalign} objects as obtained from
  calling \code{read.fasta.pdb} or two \sQuote{pdb} objects as
  obtained from calling \code{read.pdb} and does structure alignment 
  based on dynamic programming with position specific scoring matrix 
  built from residue distance. When provided with \sQuote{3dalign} 
  objects containing multiple sequences, it does profile alignment 
  in a similar way to \sQuote{CLUSTALW} program. Gap penalty is applied
  to the openning of gaps in the alignment.
}
\usage{
  straln.pair(pdbs1, pdbs2, gap.penalty = -9,
         ret.score = FALSE, fit = FALSE, model = c("ca", "backbone"),
         method = c("average", "maximal"), ncore = NULL)
}
\arguments{
  \item{pdbs1}{ A \sQuote{3dalign} object as obtained from calling \code{read.fasta.pdb}. }
  \item{pdbs2}{ A \sQuote{3dalign} object to be aligned with \sQuote{pdbs1}. }
  \item{gap.penalty}{ Penalty score for opening gaps. }
  \item{ret.score}{ logical, if TRUE return a list containing both alignment and the score. }
  \item{fit}{ logical, if TRUE doing fitting before alignment (not supported yet). }
  \item{model}{ Residue distance is calculated based on either \sQuote{CA} or \sQuote{backbone} atoms. }
  \item{method}{ Method used in profile alignment; either take average or maximal score of all sequence
                 pairs from the two alignments for each position. }
  \item{ncore}{ number of CPU cores used to do the calculation. By default, all cores will be used}
}
\value{
  if \sQuote{ret.score} is FALSE, returns a list similar to \sQuote{3dalin} class.
  if \sQuote{ret.score} is TRUE, returns:
  \item{aln}{ A list object similar to \sQuote{3dalign} class. }
  \item{score}{ The alignment score. }
}
\note{ 
  A couple things could be improved such as the support of \sQuote{fit} and separated penalty
  for gap openning and extension.
}  
\seealso{
   \code{\link{straln}}, \code{\link{seqaln}}, \code{\link{seqaln.pair}}, \code{\link{struct.aln}}
}
\examples{

\dontrun{

# Pair-wise structure alignment
#files <- get.pdb(c("1as0","5p21"), path=tempdir(), overwrite=TRUE)
#pdbs <- pdbaln(files)
#pdbs$xyz <- pdbfit(pdbs)
#pdb <- pdbs2pdb(pdbs, inds=c(1,2))
#aln <- straln.pair(pdb[[1]], pdb[[2]])

# Profile alignment

}
}
\keyword{analysis}
