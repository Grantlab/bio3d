\name{pdb.annotate}
\alias{pdb.annotate}
\title{ Get Customizable Annotations From PDB }
\description{
  Get customizable annotations for query results from PDB.  
}
\usage{
pdb.annotate(ids, anno.terms = NULL, unique = FALSE, verbose = FALSE)
}
\arguments{
  \item{ids}{ A charater vector of one or more 4-letter PDB
    codes/identifiers of the files for query. }
  \item{anno.terms}{ Terms can be used for query. The "anno.terms" can be
    "structureId", "experimentalTechnique", "resolution", "chainId",
    "ligandId", "ligandName", "source",  "scopDomain", "classification",
    "compound","title",  "citation", "citationAuthor", "journalName",
    "publicationYear".  If anno.terms=NULL, all information would be
    returned. }
  \item{unique}{ logical, if TRUE only unique PDB entries are
    returned. Alternatively data for each chain ID is provided. }
  \item{verbose}{ logical, if TRUE details of the RCurl \code{postForm}
    routine is printed. }
}
\details{
  Given a list of PDB IDs and query terms, this function will download
  the required information from PDB, and return a data frame of query
  results.
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
}
\keyword{ utilities }
