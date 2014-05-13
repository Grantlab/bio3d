\name{pdb.annotate}
\alias{pdb.annotate}
\title{ Get Customizable Annotations From PDB }
\description{
  Get customizable annotations for query results from PDB.  
}
\usage{
pdb.annotate(ids, anno.terms=NULL, unique=FALSE)
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
  \item{unique}{ logical, if TRUE the annotation is provided for unique
    PDB entries. By default the annotation is provided in the same order
    (and length) as the input list of PDB ids. }
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
\author{ Hongyang Li & Barry Grant}
\examples{

# Fetch all annotation terms
pdbanno <- pdb.annotate(c("6Q21_B", "1NVW", "1P2U_A"))

# Access terms, e.g. ligand names:
pdbanno$ligandName

# Fetch only specific terms
pdb.annotate(c("6Q21_B", "1NVW", "1P2U_A"),
             anno.terms = c("ligandId", "experimentalTechnique",
                            "publicationYear", "citation"))

}
\keyword{ utilities }