\name{pdb.annotate}
\alias{pdb.annotate}
\alias{pdb.pfam}
\title{ Get Customizable Annotations From PDB Or PFAM Databases}
\description{
  Get customizable annotations for query results from PDB or PFAM.  
}
\usage{
pdb.annotate(ids, anno.terms = NULL, unique = FALSE, verbose = FALSE)
pdb.pfam(ids, best.only = TRUE, compact = TRUE)
}
\arguments{
  \item{ids}{ A charater vector of one or more 4-letter PDB
    codes/identifiers of the files for query. }
  \item{anno.terms}{ Terms can be used for query. The "anno.terms" can be
    "structureId", "experimentalTechnique", "resolution", "chainId",
    "ligandId", "ligandName", "source",  "scopDomain", "classification",
    "compound","title",  "citation", "citationAuthor", "journalName",
    "publicationYear", "rObserved", or "rFree".  If anno.terms=NULL, 
    all information would be returned. }
  \item{unique}{ logical, if TRUE only unique PDB entries are
    returned. Alternatively data for each chain ID is provided. }
  \item{verbose}{ logical, if TRUE details of the RCurl \code{postForm}
    routine is printed. }
  \item{best.only}{ logical, if TRUE only the lowest eValue match for a 
    given input id will be reported. Otherwise all significant matches 
    will be returned. }
  \item{compact}{ logical, if TRUE only a subset of annotation terms 
    are returned. Otherwise full match details are reported (see 
    examples).}      
}
\details{
  Given a list of PDB IDs (and query terms for the pdb.annotate function),  
  these functions will download annotation information from the RCSB PDB 
  and PFAM databases.
}
\value{
  Returns a data frame of query results with a row for each PDB record,
  and annotation terms column-wise. 
}
\author{ Hongyang Li, Barry Grant, Lars Skjaerven}
\examples{
\donttest{
# PDB server connection required - testing excluded

# Fetch all annotation terms
ids <- c("6Q21_B", "1NVW", "1P2U_A")
anno <- pdb.annotate(ids)

# Access terms, e.g. ligand names:
anno$ligandName

## only unique PDB IDs
anno <- pdb.annotate(ids, unique=TRUE)

# Fetch only specific terms
pdb.annotate(ids, anno.terms = c("ligandId", "citation"))

}

\dontrun{
# PFAM server connection required - testing excluded

# Find PFAM annotations of PDB entries
pdb.pfam("6Q21_A", "1NVW", "1P2U_A")

# More details and a not fond entry warning
pdb.pfam(c("1P2U_A", "6Q21_B"), compact=FALSE)
}

}
\keyword{ utilities }
