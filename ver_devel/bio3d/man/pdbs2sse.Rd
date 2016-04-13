\name{pdbs2sse}
\alias{pdbs2sse}
\title{ SSE annotation for a PDBs Object }
\description{
  Returns secondary structure element (SSE) annotation (\code{"sse"} object) for a structure in the provided \code{"pdbs"} object.
}
\usage{
pdbs2sse(pdbs, ind = 1, rm.gaps = TRUE, resno = TRUE, pdb = FALSE, \dots)
}
\arguments{
  \item{pdbs}{ a list of class \code{"pdbs"} containing PDB file
    data, as obtained from \code{read.fasta.pdb} or
    \code{pdbaln}. }
  \item{ind}{ numeric index pointing to the PDB in which the SSE should
    be provided. }
  \item{rm.gaps}{ logical, if TRUE SSEs spanning gap containing columns are
    omitted from the output in the resulting \code{sse} object. }
  \item{resno}{ logical, if TRUE output is in terms of residue numbers
    rather than residue index (position in sequence). }
  \item{pdb}{ logical, if TRUE function \code{dssp} will be called on
    the corresponding \code{pdb} object rather than to use \code{pdbs$sse} to
    obtain the SSE object. }
   \item{\dots}{ arguments passed to function \code{dssp}. }
}
\details{
  This function provides a \code{"sse"} list object containing
  secondary structure elements (SSE) annotation data for a particular
  structure in the provided \code{"pdbs"} object. Residue numbers are
  provided relative to the alignment in the \code{"pdbs"} object.

  See examples for more details.
}
\value{
  Returns a list object of class \code{sse}.
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Lars Skjaerven }
\seealso{ \code{\link{dssp}}, \code{\link{pdbaln}},
  \code{\link{read.fasta.pdb}}. }
\examples{
\dontrun{
attach(transducin)

## calculate RMSF
rf <- rmsf(pdbs$xyz)

## Fetch SSE annotation, output in terms of alignment index
sse <- pdbs2sse(pdbs, ind=1, rm.gaps=FALSE, resno=FALSE)

## Add SSE annotation to plot
plotb3(rf, sse=sse)


## Calculate RMSF only for non-gap columns
gaps.pos <- gap.inspect(pdbs$xyz)
rf <- rmsf(pdbs$xyz[, gaps.pos$f.inds])

## With gap columns removed, output in terms of residue number
sse <- pdbs2sse(pdbs, ind=1, rm.gaps=TRUE, resno=TRUE)
gaps.res <- gap.inspect(pdbs$ali)
plotb3(rf, sse=sse, resno=pdbs$resno[1, gaps.res$f.inds])

}
}
\keyword{ utilities }
