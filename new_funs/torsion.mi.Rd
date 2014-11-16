\name{torsion.mi}
\alias{torsion.mi}
\title{ Mutual Information Cross-correlation of Sidechain Torsion Angles }
\description{
  Calculate mutual information cross-correlation based on sidechain 
  torsion angle chi1.
}
\usage{
torsion.mi(xyz, pdb, normalize = TRUE, 
    xyz2tor = TRUE, ncore = NULL)
}
\arguments{
  \item{xyz}{ Cartesian coordinates/torsion angles as a MxN matrix. }
  \item{pdb}{ A PDB structure object as obtained from \code{\link{read.pdb}}. }
  \item{normalize}{logical, if TRUE normalize return values to the 
      range [0,1]. }
  \item{xyz2tor}{logical, if TRUE call \code{\link{xyz2torsion}} to convert 
      \code{xyz} to torsion angles. }
  \item{ncore}{ Number of CPU cores used to do the calculation. By default (NULL), use
    all detected CPU cores. }
}
\details{
  The Cartesian coordinate matrix is first converted to sidechain torsion angles (chi1)
  for each residue (neglecting alanine and glycine) with function \code{\link{xyz2torsion}} (with \code{xyz2tor} = TRUE).
  Three states are defined for each type of amino acid with valid chi1 except for proline,
  with chi1 values being:
  \itemize{
     \item gauche+, [0, 120) degree.
     \item gauche-, [-120, 0) degree.
     \item trans, others.
  }
  Two states are defined for proline: gauche+, [0, 90) or [-180, -90), and gauche-, others.

  Mutual information for each pair of residues is then calculated with the formula,
     \deqn{I(X, Y) = H(X) + H(Y) - H(X, Y),}
  where \deqn{H(X) = - \sum\limit_i P(X=i)*log(P(X=i))} is the Shannon entropy for random variable X.
  Cross-correlation is obtained by normalizing I(X, Y), 
     \deqn{C(X, Y) = 2*I(X, Y)/(H(X) + H(Y)),}
  with the fact that
     \deqn{I(X, Y) <= (H(X) + H(Y)) / 2.}
}
\value{
  Returns a matrix of residue-residue mutual information or normalized cross-correlation. 
}
\note{
  Candidates for the denominator to normalize I(X, Y) include
     \deqn{
     D1 = min{H(X), H(Y)}, 
     D2 = (H(X), H(Y))/2, 
     D3 = max{H(X), H(Y)}. 
     }
  We choose D2, but other possibility can not be excluded.
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
  Bowman, G.R. and Geissler, P.L. (2012) \emph{PNAS} \bold{109}, 11681-11686.
}
\author{ Xin-Qiu Yao, Guido Scarabelli }
\seealso{ \code{\link{torsion.xyz}}, \code{\link{torsion.pdb}},
  \code{\link{xyz2torsion}}. }
\examples{
\dontrun{
pdb <- read.pdb("sys_gdp_nowat.pdb")
xyz <- read.ncdf("dyna.01.traj.nc")
cij <- torsion.mi(xyz, pdb)

# plot DCCM
plot(cij)
}
}
\keyword{ utilities }
