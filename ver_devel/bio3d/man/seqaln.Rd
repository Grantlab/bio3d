\name{seqaln}
\alias{seqaln}
\title{ Sequence Alignment with MUSCLE}
\description{
  Create multiple alignments of amino acid or nucleotide sequences
  according to the method of Edgar.
}
\usage{
seqaln(aln, id=NULL, profile=NULL, exefile="muscle", outfile="aln.fa", 
       protein=TRUE, seqgroup=FALSE, refine=FALSE, extra.args="",
       verbose=FALSE)  
}
\arguments{
  \item{aln}{ a sequence character matrix, as obtained from
    \code{\link{seqbind}}, or an alignment list object as obtained from
      \code{\link{read.fasta}}. }
  \item{id}{ a vector of sequence names to serve as sequence
    identifers. }
  \item{profile}{ a profile alignment of class \sQuote{fasta}
    (e.g. obtained from \code{\link{read.fasta}}). The alignment
    \code{aln} will be added to the profile. }
  \item{exefile}{ file path to the \sQuote{MUSCLE} program on your system (i.e.
    how is \sQuote{MUSCLE} invoked). Alternatively, \sQuote{CLUSTALO} can
    be used. }
  \item{outfile}{ name of \sQuote{FASTA} output file to which alignment
    should be written. }
  \item{protein}{ logical, if TRUE the input sequences are assumed to be
    protein not DNA or RNA. }
  \item{seqgroup}{ logical, if TRUE similar sequences are grouped
    together in the output. }
  \item{refine}{ logical, if TRUE the input sequences are assumed to already
  be aligned, and only tree dependent refinement is performed. }
  \item{extra.args}{ a single character string containing extra command
    line arguments for the alignment program. }
  \item{verbose}{ logical, if TRUE \sQuote{MUSCLE} warning and error
    messages are printed. } 
}
\details{
  Sequence alignment attempts to arrange the sequences of protein,
  DNA or RNA, to highlight regions of shared similarity that may
  reflect functional, structural, and/or evolutionary relationships
  between the sequences.

  Aligned sequences are represented as rows within a matrix. Gaps
  (\sQuote{-}) are inserted between the aminoacids or nucleotides so that
  equivalent characters are positioned in the same column.

  This function calls the \sQuote{MUSCLE} program, to perform a multiple
  sequence alignment, which MUST BE INSTALLED on your system and in the
  search path for executables.
 
  If you have a large number of input sequences (a few thousand), or they are
  very long, the default settings may be too slow for practical
  use. A good compromise between speed and accuracy is to run just the
  first two iterations of the \sQuote{MUSCLE} algorithm by setting the
  \code{extra.args} argument to \dQuote{-maxiters 2}.

  You can set \sQuote{MUSCLE} to improve an existing alignment by setting
  \code{refine} to TRUE.

  To inspect the sequence clustering used by \sQuote{MUSCLE} to produce
  alignments, include \dQuote{-tree2 tree.out} in the \code{extra.args}
  argument. You can then load the \dQuote{tree.out} file with the
  \sQuote{read.tree} function from the \sQuote{ape} package.

  \sQuote{CLUSTALO} can be used as an alternative to \sQuote{MUSCLE} by
  specifiying \code{exefile='clustalo'}. This might be useful e.g. when
  adding several sequences to a profile alignment. 
}
\value{
  A list with two components:
  \item{ali}{ an alignment character matrix with a row per sequence and
    a column per equivalent aminoacid/nucleotide. }
  \item{ids}{ sequence names as identifers.}
  \item{call}{ the matched call. }
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
  
  \sQuote{MUSCLE} is the work of Edgar:
  Edgar (2004) \emph{Nuc. Acid. Res.} \bold{32}, 1792--1797.
  
  Full details of the \sQuote{MUSCLE} algorithm, along with download and
  installation instructions can be obtained from:\cr
  \url{http://www.drive5.com/muscle}.
}
\author{ Barry Grant }
\note{
  A system call is made to the \sQuote{MUSCLE} program, which must be
  installed on your system and in the search path for executables.
}
\seealso{ \code{\link{read.fasta}}, \code{\link{read.fasta.pdb}},
  \code{\link{get.seq}}, \code{\link{seqbind}}, 
  \code{\link{pdbaln}}, \code{\link{plot.fasta}},
  \code{\link{blast.pdb}}
}
\examples{
\dontrun{
##-- Basic sequence alignemnt
seqs <- get.seq(c("4q21_A", "1ftn_A"))
aln <- seqaln(seqs)

##-- add a sequence to the (profile) alignment
seq <- get.seq("1tnd_A")
aln <- seqaln(seq, profile=aln)


##-- Read a folder/directory of PDB files
#pdb.path <- "my_dir_of_pdbs"
#files  <- list.files(path=pdb.path ,
#                     pattern=".pdb",
#                     full.names=TRUE)

##-- Use online files
files <- get.pdb(c("4q21","1ftn"), URLonly=TRUE)

##-- Extract and store sequences
raw <- NULL
for(i in 1:length(files)) {
  pdb <- read.pdb(files[i])
  raw <- seqbind(raw, pdbseq(pdb) )
}

##-- Align these sequences
aln <- seqaln(raw, id=files, outfile="seqaln.fa")

##-- Read Aligned PDBs storing coordinate data
pdbs <- read.fasta.pdb(aln) 

## Sequence identity
seqidentity(aln)

## Note that all the above can be done with the pdbaln() function:
#pdbs <- pdbaln(files)


##-- For identical sequences with masking use a custom matrix
aln <- list(ali=seqbind(c("X","C","X","X","A","G","K"),
                        c("C","-","A","X","G","X","X","K")),
            id=c("a","b"))

temp <- seqaln(aln=aln, outfile="temp.fas", protein=TRUE,
               extra.args= paste("-matrix",
                system.file("matrices/custom.mat", package="bio3d"),
                "-gapopen -3.0 ",
                "-gapextend -0.5",
                "-center 0.0") )

}
}
\keyword{ utilities }
