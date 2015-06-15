\name{pfam}
\alias{pfam}
\title{ Download Pfam FASTA Sequence Alignment }
\description{
  Downloads FASTA sequence alignment from the Pfam database.
}
\usage{
pfam(id, alignment = "seed", verbose = FALSE)
}
\arguments{
  \item{id}{ the Pfam familiy identifier (e.g \sQuote{Piwi}) or accession
    (e.g. \sQuote{PF02171}). } 
  \item{alignment}{ the alignment type. Allowed values are:
    \sQuote{seed}, \sQuote{ncbi}, \sQuote{full},
    \sQuote{metagenomics}. }
  \item{verbose}{ logical, if TRUE details of the download process
    is printed. }
}
\details{
  This is a basic function to download a multiple sequence alignment for
  a protein family from the Pfam database. 
}
\value{
  A \sQuote{fasta} object with the following components:
  \item{ali }{ an alignment character matrix with a row per sequence and
    a column per equivalent aminoacid/nucleotide. }
  \item{ids }{ sequence names as identifiers. }
  \item{call }{ the matched call. }
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Lars Skjaerven }
\note{
  Full more information on the Pfam database:\cr
  \url{http://pfam.xfam.org}
}
\seealso{ \code{\link{read.fasta}},
  \code{\link{hmmer}}, \code{\link{get.seq}},
  \code{\link{uniprot}}
}
\examples{
\donttest{
# PFAM server connection required - testing excluded

aln <- pfam("piwi")
aln <- pfam("PF02171")

seq <- get.seq("1rx2_A", outfile = tempfile())
hmm <- hmmer(seq, type="hmmscan", db="pfam")
aln <- pfam(hmm$acc[1])

# Or much more simply for RCSB PDB entries:
acc <- pdb.pfam("1rx2_A", compact=FALSE)$pfamAcc
aln <- pfam(acc)
}
}
\keyword{ utilities }
