\name{sse.vector}
\alias{sse.vector}
\title{
    Output Secondary Structure Vectors
}
\description{
    Output a user defined secondary structure vector for a given PDB input structure.
}
\usage{
    sse.vector(pdb, sse = pdb, coil = "-", helix = "H", sheet = "E", calpha=TRUE)
}
\arguments{
  \item{pdb}{ a PDB structure object obtained from \code{\link{read.pdb}}. }
  \item{sse}{ a secondary structure object as returned from \code{\link{read.pdb}},
    \code{\link{dssp}} or \code{\link{stride}}. }
  \item{coil}{ a single element vector for representing non helix and sheet regions. }
  \item{helix}{ a single element vector for representing alpha helices. }
  \item{sheet}{ a single element vector for representing beta strands. }
  \item{sheet}{ logical, if TRUE a per-residue (calpha atom) vector will be returned. If FALSE a per atom vector will be returned. See examples section for details. }  
}
\details{
    Can be used to output a per-residue secondary structure sequence vector or a per-atom color vector for visualization. See examples section.
}
\value{
Returns a vector of user defined coil, helix and sheet values per residue or per atom of the input PDB structure. The names attribute of this returned vector corresponds to residue number and chain identifiers for each position or atom.
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Barry Grant }
\seealso{ \code{\link{read.pdb}}, \code{\link{dssp}}, \code{\link{stride}}, \code{\link{view.pdb}} }
\examples{

\dontrun{

pdb <- read.pdb("4q21")
## Per-residue secondary structure vector
paste( sse.vector(pdb), collapse="")

## Per-atom color vector for secondary structure visualization
myCols <- sse.vector(pdb, coil="gray", helix="green", sheet="cyan", calpha=FALSE)
view(pdb, col=myCols, mode="calpha")

## Compare PDB, DSSP and STRIDE defined secondary structure
cbind( sse.vector(pdb), sse.vector(pdb, dssp(pdb)), sse.vector(pdb, stride(pdb)) )

}

}
\keyword{ utilities }
