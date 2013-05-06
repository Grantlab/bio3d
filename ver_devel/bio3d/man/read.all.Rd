\name{read.all}
\alias{read.all}
\title{ Read Aligned Structure Data}
\description{
  Read aligned PDB structures and store their equalvalent atom data, including xyz
  coordinates, residue numbers, residue type and B-factors.
}
\usage{
read.all(aln, pdb.path = "", pdbext = "", sel = NULL, ...)
}

\arguments{
  \item{aln}{ an alignment data structure obtained with
    \code{\link{read.fasta}}. }
  \item{pdb.path}{ path to PDB files. }
  \item{pdbext}{ the file name extention of the PDB files. }
  \item{sel}{ a selection string detailing the atom type data to store
    (see function store.atom) }
  \item{\dots}{ other parameters for \code{\link{read.pdb}}. }   
}
\details{
  The input \code{aln}, produced with \code{\link{read.fasta}}, must
  have identifers (i.e. sequence names) that match the PDB file
  names. For example the sequence corresponding to the structure
  file \dQuote{mypdbdir/1bg2.pdb} should have the identifer
  \sQuote{mypdbdir/1bg2.pdb} or \sQuote{1bg2} if input \sQuote{pdb.path}
  and \sQuote{pdbext} equal \sQuote{mypdbdir/} and \sQuote{pdb}. See the
  examples below.

  Sequence miss-matches will generate errors.  Thus, care should be taken
  to ensure that the sequences in the alignment match the sequences in
  their associated PDB files.
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
  \item{all }{numeric matrix of aligned equalvelent atom coordinates. }
  \item{all.elety}{numeric matrix of aligned atom element types. }
  \item{all.resid}{numeric matrix of aligned three-letter residue codes. }
  \item{all.resno}{numeric matrix of aligned residue numbers. }
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Barry Grant }
\note{
  This function is still in development and is NOT part of the
  offical bio3d package.
  
  The sequence character \sQuote{X} is useful for masking unusual
  or unknown residues, as it can match any other residue type.
}
\seealso{
  \code{\link{read.fasta}}, \code{\link{read.pdb}},
  \code{\link{core.find}}, \code{\link{fit.xyz}}  }
\examples{
# still working on speeding this guy up
\dontrun{
# Read sequence alignment
file <- system.file("examples/kif1a.fa",package="bio3d")
aln  <- read.fasta(file)

# Read aligned PDBs
pdb.path = system.file("examples",package="bio3d")
sel <- c("N", "CA", "C", "O", "CB", "*G", "*D",  "*E", "*Z")
pdbs <- read.all(aln, pdb.path = pdb.path, pdbext = ".ent", sel=sel)

atm <- colnames(pdbs$all)
ca.ind  <- which(atm == "CA")
core <- core.find(pdbs)
core.ind <- c( matrix(ca.ind, nrow=3)[,core$c0.5A.atom] )

nxyz <- fit.xyz(pdbs$all[1,], pdbs$all,
               fixed.inds  = core.ind,
               mobile.inds = core.ind)

ngap.col <- gap.inspect(nxyz)
npc.xray <- pca.xyz(nxyz[ ,ngap.col$f.inds])


a <- mktrj.pca(npc.xray, pc=1, file="pc1-all.pdb",
               elety=pdbs$all.elety[1,unique( ceiling(ngap.col$f.inds/3) )],
               resid=pdbs$all.resid[1,unique( ceiling(ngap.col$f.inds/3) )],
               resno=pdbs$all.resno[1,unique( ceiling(ngap.col$f.inds/3) )] )


}

}
\keyword{ IO }
