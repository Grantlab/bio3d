\name{check.utility}
\alias{check.utility}
\title{ Check on Missing Utility Programs }
\description{
  Internally used in examples, tests, or vignettes.
}
\usage{
  check.utility(x = c("muscle", "dssp", "stride", "mustang", "makeup"), quiet = TRUE)
}
\arguments{
  \item{x}{ Names of one or more utility programs to check. }
  \item{quiet}{ logical, if TRUE no warning or message printed. }
}
\details{
  Check if requested utility programs are availabe or not.
}
\value{
  logical, TRUE if programs are available and FALSE if any one of them is missing.
}
\examples{
   check.utility(c("muscle", "dssp"), quiet=FALSE)
   if(!check.utility("mustang")) 
      cat(" The utility program, MUSTANG, is missing on your system\n")
}
\keyword{ utilities }
