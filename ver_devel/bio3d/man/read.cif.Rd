\name{read.cif}
\alias{read.cif}
\title{ Read mmCIF File}
\description{
  Read a Protein Data Bank (mmCIF) coordinate file.
}
\usage{
read.cif(file, maxlines = -1, multi = FALSE,
         rm.insert = FALSE, rm.alt = TRUE, verbose = TRUE)
}
\arguments{
  \item{file}{ a single element character vector containing the name of
    the mmCIF file to be read, or the four letter PDB identifier for
    online file access. }
  \item{maxlines}{ the maximum number of lines to read before giving up
    with large files. By default if will read up to the end of input on
    the connection. }
  \item{multi}{ logical, if TRUE multiple ATOM records are read for all
    models in multi-model files and their coordinates returned. }
  \item{rm.insert}{ logical, if TRUE PDB insert records are ignored. }
  \item{rm.alt}{ logical, if TRUE PDB alternate records are ignored. }
  \item{verbose}{ print details of the reading process. }
}
\details{
  The current version of \code{read.cif} reads only ATOM/HETATM records
  and creates a \code{pdb} object of the data.

  See \code{read.pdb} for more info.
}
\value{
  Returns a list of class \code{"pdb"} with the following components:
  \item{atom}{ a data.frame containing all atomic coordinate ATOM and 
    HETATM data, with a row per ATOM/HETATM and a column per record type.  
    See below for details of the record type naming convention (useful for 
    accessing columns). }
  \item{xyz }{ a numeric matrix of class \code{"xyz"} containing the
    ATOM and HETATM coordinate data. }
  \item{calpha }{ logical vector with length equal to \code{nrow(atom)}
    with TRUE values indicating a C-alpha \dQuote{elety}. }
  \item{call }{ the matched call. }
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Lars Skjaerven }
\seealso{
  \code{\link{read.pdb}}
  \code{\link{atom.select}}, \code{\link{write.pdb}},
  \code{\link{trim.pdb}}, \code{\link{cat.pdb}}, 
  \code{\link{read.prmtop}}, \code{\link{as.pdb}},
  \code{\link{read.dcd}}, \code{\link{read.ncdf}}, 
}
\examples{
## Read a mmCIF file from the RCSB online database
# cif <- read.cif("1hel")
}
\keyword{ IO }

