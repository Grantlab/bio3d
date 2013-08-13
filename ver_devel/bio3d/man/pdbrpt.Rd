\name{pdbrpt}
\alias{pdbrpt}
\title{ Get Customizable Reports From PDB }
\description{
  Get customizable reports for query results from PDB.  
}
\usage{
pdbrpt(ids, anno.terms=NULL)
}
\arguments{
 \item{ids}{ A charater vector of one or more 4-letter PDB codes/identifiers of the files for query. }
 \item{anno.terms}{ Terms can be used for query. The "anno.terms" can be "structureId", "experimentalTechnique", "resolution", "chainId", "ligandId", "ligandName", "source", "scopDomain", "classification", "compound","title", "citation", "citationAuthor", "journalName", "publicationYear". If anno.terms=NULL, all information would be returned.
}
}
\details{
 Given a list of PDB IDs and query terms, this function will download the required information from PDB, remove redundant information and return a matrix of query results. Require XML package. 
}
\value{
 Return a matrix of query results.
}
\author{ Hongyang Li & Barry Grant}
\examples{
pdbrpt(c("6Q21_B", "1NVW", "1P2U_A"), anno.terms = c("ligandI", "exp", "pub", "citation"))
## Return a matrix of the following information:
# PDB ID: 6Q21, 1NVW, 1P2U
# Fields: ligandId, experimentalTechnique, publicationYear, citation
}
