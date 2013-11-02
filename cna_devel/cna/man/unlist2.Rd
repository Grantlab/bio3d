\name{unlist2}

\alias{unlist2}

\title{A replacement for unlist() that does not mangle the names}

\description{
  \code{unlist2} is a replacement for \code{base::unlist()} that
  does not mangle the names.
}

\usage{
  unlist2(x, recursive=TRUE, use.names=TRUE, what.names="inherited")
}

\arguments{
  \item{x, recursive, use.names}{See \code{?unlist}.}
  \item{what.names}{\code{"inherited"} or \code{"full"}.}
}

\details{
  Use this function if you don't like the mangled names returned
  by the standard \code{unlist} function from the base package.
  Using \code{unlist} with annotation data is dangerous and it is
  highly recommended to use \code{unlist2} instead.
}

\author{Herve Pages}

\seealso{
  \code{\link{unlist}}
}

\examples{
  x <- list(A=c(b=-4, 2, b=7), B=3:-1, c(a=1, a=-2), C=list(c(2:-1, d=55), e=99))
  unlist(x)
  unlist2(x)

  library(hgu95av2.db)
  egids <- c("10", "100", "1000")
  egids2pbids <- mget(egids, revmap(hgu95av2ENTREZID))
  egids2pbids

  unlist(egids2pbids)   # 1001, 1002, 10001 and 10002 are not real
                        # Entrez ids but are the result of unlist()
                        # mangling the names!

  unlist2(egids2pbids)  # much cleaner! yes the names are not unique
                        # but at least they are correct...
}

\keyword{utilities}
\keyword{manip}
