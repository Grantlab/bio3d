\name{write.ncdf}
\alias{write.ncdf}
\title{ Write AMBER Binary netCDF files }
\description{
  Write coordinate data to a binary netCDF trajectory file.
}
\usage{
write.ncdf(x, trjfile = "R.ncdf")
}
\arguments{
  \item{x}{ A numeric matrix of xyz coordinates with a frame/structure
    per row and a Cartesian coordinate per column. }
  \item{trjfile}{ name of the output trajectory file. }
}
\details{
  Writes an AMBER netCDF (Network Common Data Form) format trajectory
  file with the help of David W. Pierce's (UCSD) ncdf package available
  from CRAN.
}
\value{
  Called for its effect.
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
  \url{http://www.unidata.ucar.edu/packages/netcdf/}
  \url{http://cirrus.ucsd.edu/~pierce/ncdf/}
  \url{http://amber.scripps.edu/netcdf/nctraj.html}
}
\author{ Barry Grant }
\note{ 
  See AMBER documentation for netCDF format description.

  NetCDF binary trajectory files are supported by the AMBER modules
  sander, pmemd and ptraj. Compared to formatted trajectory files, the
  binary trajectory files are smaller, higher precision and
  significantly faster to read and write.

  NetCDF provides for file portability across architectures, allows for
  backwards compatible extensibility of the format and enables the files
  to be self-describing. Support for this format is available in VMD.
}
\seealso{ \code{\link{read.dcd}}, \code{\link{read.ncdf}},
  \code{\link{read.pdb}}, \code{\link{write.pdb}},
  \code{\link{atom.select}} }
\examples{
\dontrun{
##-- Read example trajectory file
trtfile <- system.file("examples/hivp.dcd", package="bio3d")
trj <- read.dcd(trtfile)

## Write to netCDF format
write.ncdf(trj, "newtrj.nc")

## Read trj
trj <- read.ncdf("newtrj.nc")
}
}
\keyword{ IO }
