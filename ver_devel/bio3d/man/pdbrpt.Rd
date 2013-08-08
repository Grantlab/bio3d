\name{pdbrpt}
\alias{pdbrpt}
\title{ Get Customizable Reports From PDB }
\description{
  Get customizable reports for query results from PDB in csv format.  
}
\usage{
pdbrpt(ids, file="pdbrpt.csv", terms=NULL)
}
\arguments{
 \item{ids}{ A charater vector of one or more 4-letter PDB codes/identifiers of the files for query. }
 \item{file}{ Output report file. }
 \item{terms}{ Terms can be used for query. Terms can be "structureId", "experimentalTechnique", "resolution", "chainId", "ligandId", "ligandName", "source", "scopDomain", "classification", "compound","title", "citation", "citationAuthor", "journalName", "publicationYear"
}
}
\details{
 Given a list of PDB IDs and query terms, this function will download the required information from PDB, remove redundant information and generate a report in csv format. Require XML package. 
}
\value{
 No returns
}
\author{ Hongyang Li & Barry Grant}
\examples{
pdbrpt(c("6Q21_B", "1NVW", "1P2U_A"), file="pdbrpt.csv", terms = c("ligandI", "exp", "pub", "citation"))
## Generate a report "pdbrpt.csv" with the following information:
# PDB ID: 6Q21, 1NVW, 1P2U
# Fields: ligandId, experimentalTechnique, publicationYear, citation
}
