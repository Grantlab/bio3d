\name{core.cmap}
\alias{core.cmap}
\title{ Identification of Contact Map Core Positions }
\description{
  Find core positions that have the largest number
  of contact with neighboring residues.
}
\usage{
core.cmap(pdbs, write.pdb = FALSE, outfile="core.pdb",
          cutoff = NULL, refine = FALSE, ncore = NULL, \dots)
}
\arguments{
  \item{pdbs}{ an alignment data structure of class \sQuote{pdbs} 
    as obtained with \code{\link{read.fasta.pdb}} or \code{\link{pdbaln}}, 
    or a numeric matrix of aligned C-alpha xyz Cartesian coordinates. }
  \item{write.pdb}{ logical, if TRUE core coordinate files, containing
    only core positions for each iteration, are written to a location
    specified by \code{outpath}. }
  \item{outfile}{ character string specifying the output directory when
    \code{write.pdb} is \sQuote{TRUE}. }  
  \item{cutoff}{ numeric value speciyfing the inclusion criteria for
    core positions. }
  \item{refine}{ logical, if TRUE explore core positions determined by
    multiple eigenvectors. By default only the eigenvector describing
    the largest variation is used. }
  \item{ncore }{ number of CPU cores used to do the calculation.
    By default (\code{ncore=NULL}) use all cores detected. }
  \item{\dots}{ arguments passed to and from functions. }
}
\details{
  This function calculates eigenvector centrality of the weighted
  contact network built based on input structure data and uses it
  to determine the core positions. 
}
\value{
  Returns a list of class \code{"select"} containing \sQuote{atom} and
  \sQuote{xyz} indices.
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Xin-Qiu Yao }
\seealso{ \code{\link{read.fasta.pdb}}, \code{\link{core.find}},
  \code{\link{fit.xyz}} }
\examples{
\dontrun{
##-- Generate a small kinesin alignment and read corresponding structures
pdbfiles <- get.pdb(c("1bg2","2ncd","1i6i","1i5s"), URLonly=TRUE)
pdbs <- pdbaln(pdbfiles)

##-- Find 'core' positions
core <- core.cmap(pdbs)
xyz <- pdbfit(pdbs, core, outpath="corefit_structures")
}
}
\keyword{ utilities }
