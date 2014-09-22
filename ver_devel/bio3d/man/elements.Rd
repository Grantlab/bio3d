\name{elements}
\alias{elements}
\docType{data}
\title{Periodic Table of the Elements}
\description{This data set gives various information on chemical elements.}
\usage{elements}
\format{
  A data frame containing for each chemical element the following information.
  \describe{
    \item{\code{num}}{atomic number}
    \item{\code{symb}}{elemental symbol}
    \item{\code{areneg}}{Allred and Rochow electronegativity  (0.0 if unknown)}
    \item{\code{rcov}}{covalent radii (in Angstrom) (1.6 if unknown)}
    \item{\code{rbo}}{"bond order" radii}
    \item{\code{rvdw}}{van der Waals radii (in Angstrom) (2.0 if unknown)}
    \item{\code{maxbnd}}{maximum bond valence (6 if unknown)}
    \item{\code{mass}}{IUPAC recommended atomic masses (in amu)}
    \item{\code{elneg}}{Pauling electronegativity (0.0 if unknown)}
    \item{\code{ionization}}{ionization potential (in eV) (0.0 if unknown)}
    \item{\code{elaffinity}}{electron affinity (in eV) (0.0 if unknown)}
    \item{\code{red}}{red value for visualization}
    \item{\code{green}}{green value for visualization}
    \item{\code{blue}}{blue value for visualization}
    \item{\code{name}}{element name}
  }
}

\source{
Open Babel (2.3.1) file: element.txt\cr

  Created from the Blue Obelisk Cheminformatics Data Repository\cr
  Direct Source: http://www.blueobelisk.org/\cr
  http://www.blueobelisk.org/repos/blueobelisk/elements.xml includes furhter bibliographic citation information\cr

  - Allred and Rochow Electronegativity from http://www.hull.ac.uk/chemistry/electroneg.php?type=Allred-Rochow\cr
  - Covalent radii from http://dx.doi.org/10.1039/b801115j\cr
  - Van der Waals radii from http://dx.doi.org/10.1021/jp8111556\cr
}

\examples{
data(elements)
elements

# Get the mass of some elements
symb <- c("C","O","H")
elements[match(symb,elements[,"symb"]),"mass"]

# Get the van der Waals radii of some elements
symb <- c("C","O","H")
elements[match(symb,elements[,"symb"]),"rvdw"]

}
\keyword{datasets}
