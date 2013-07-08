\name{read.ncdf}
\alias{read.ncdf}
\title{ Read AMBER Binary netCDF files }
\description{
  Read coordinate data from a binary netCDF trajectory file.
}
\usage{
read.ncdf(trjfile, headonly = FALSE, verbose = TRUE, time = FALSE,
          first = NULL, last = NULL, stride = 1, cell = FALSE)
}
\arguments{
  \item{trjfile}{ name of trajectory file to read. A vector if treat
    a batch of files }
  \item{headonly}{ logical, if TRUE only trajectory header information
    is returned. If FALSE only trajectory coordinate data is returned. }
  \item{verbose}{ logical, if TRUE print details of the reading process. }
  \item{time}{ logical, if TRUE the \code{first} and \code{last} have the
    time unit ps; Otherwise the unit is the frame number. }
  \item{first}{ starting time or frame number to read; If NULL, start
    from the begining of the file(s). }
  \item{last}{ read data until \code{last} time or frame number; If NULL or
    equal to -1, read until the end of the file(s). }
  \item{stride}{ take at every \code{stride} frame(s) }
  \item{cell}{ logical, if TRUE and \code{headonly} is FALSE return cell 
    information only. Otherwise, return header or coordinates.}
}
\details{
  Reads a AMBER netCDF format trajectory file with the help of David
  W. Pierce's (UCSD) ncdf package available from CRAN.
}
\value{
   A list of trajectory header data,  a numeric matrix of xyz
   coordinates with a frame/structure per row and a Cartesian 
   coordinate per column, or a numeric matrix of cell information 
   with a frame/structure per row and lengths and angles per column.
   If time=TRUE, row names of returned coordinates or cell are
   set to be the physical time of corresponding frames.
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
  \url{http://www.unidata.ucar.edu/packages/netcdf/}
  \url{http://cirrus.ucsd.edu/~pierce/ncdf/}
  \url{http://ambermd.org/netcdf/nctraj.html}
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

  If you experience problems reading your trajectory file with read.ncdf()
  consider first reading your file into VMD and from there exporting a
  new DCD trajectory file with the 'save coordinates' option. This new
  file should be easily read with read.dcd().
}
\seealso{ \code{\link{read.dcd}}, \code{\link{write.ncdf}},
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
