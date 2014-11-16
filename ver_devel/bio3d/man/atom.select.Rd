\name{atom.select}
\alias{atom.select}
\alias{print.select}
\title{ Atom Selection }
\description{
  Return the \sQuote{atom} and \sQuote{xyz} coordinate indices of a \sQuote{pdb} structure
  object corresponding to the intersection of a hierarchical selection.
}
\usage{
atom.select(\dots)
\method{print}{select}(x, \dots)
}
\arguments{
  \item{\dots}{ arguments passed to \code{atom.select.pdb},
    \code{atom.select.prmtop}, or \code{print}. }
  \item{x}{ a atom.select object as obtained from \code{\link{atom.select}}.  }
}
\details{
  See \code{atom.select.pdb} or \code{atom.select.prmtop} for more
  information. 
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Barry Grant }
\keyword{utilities}

