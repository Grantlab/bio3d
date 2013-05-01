\name{core.find}
\alias{core.find}
\title{ Identification of Invariant Core Positions }
\description{
  Perform iterated rounds of structural superposition to identify the
  most invariant region in an aligned set of protein structures.
}
\usage{
core.find(aln, shortcut = FALSE, rm.island = FALSE,
          verbose = TRUE, stop.at = 15, stop.vol = 0.5,
          write.pdbs = FALSE, outpath="core_pruned/",
          ncore = 1, nseg.scale = 1)
}
\arguments{
  \item{aln}{ a numeric matrix of aligned C-alpha xyz Cartesian
    coordinates. For example an alignment data structure obtained with
    \code{\link{read.fasta.pdb}} or a trajectory subset obtained from \code{\link{read.dcd}}. }
  \item{shortcut}{ if TRUE, remove more than one position at a time. }
  \item{rm.island}{ remove isolated fragments of less than three
    residues. }
  \item{verbose}{ logical, if TRUE a \dQuote{core\_pruned} directory
    containing \sQuote{core structures} for each iteraction is written 
    to the current directory. }
  \item{stop.at}{ minimal core size at which iterations should be
    stopped. }
  \item{stop.vol}{ minimal core volume at which iterations should be
    stopped. }  
  \item{write.pdbs}{ logical, if TRUE core coordinate files, containing
    only core positions for each iteration, are written to a location
    specified by \code{outpath}. }
  \item{outpath}{ character string specifying the output directory when
    \code{write.pdbs} is TRUE. }  
  \item{ncore }{ number of CPU cores used to do the calculation.
    \code{ncore>1} requires multicore package installed. }
  \item{nseg.scale }{ split input data into specified number of segments
    prior to running multiple core calculation. See \code{\link{fit.xyz}}.
  }
}
\details{
  This function attempts to iteratively refine an initial structural
  superposition determined from a multiple alignment.
  This involves iterated rounds of superposition, where at each round the
  position(s) displaying the largest differences is(are) excluded from the
  dataset. 
  The spatial variation at each aligned position is determined from the
  eigenvalues of their Cartesian coordinates (i.e. the variance of the
  distribution along its three principal directions). Inspired by the
  work of Gerstein \emph{et al.} (1991, 1995), an ellipsoid of
  variance is determined from the eigenvalues, and its volume is taken as
  a measure of structural variation at a given position.

  Optional \dQuote{core PDB files} containing core positions, upon which
  superposition is based, can be written to a location specified by
  \code{outpath} by setting \code{write.pdbs=TRUE}.  These files are
  useful for examining the core filtering process by visualising them in a
  graphics program.
  
 }
\value{
  Returns a list of class \code{"core"} with the following components:
  \item{volume }{total core volume at each fitting iteration/round.}
  \item{length }{core length at each round.}
  \item{resno }{residue number of core residues at each round (taken
    from the first aligned structure) or, alternatively, the numeric
    index of core residues at each round.}
  \item{atom }{atom indices of core atoms at each round.}
  \item{xyz }{xyz indices of core atoms at each round.}
  \item{c1A.atom }{atom indices of core positions with a total volume
    under 1 Angstrom\^3.}
  \item{c1A.xyz }{xyz indices of core positions with a total volume
    under 1 Angstrom\^3.}
  \item{c1A.resno }{residue numbers of core positions with a total volume
    under 1 Angstrom\^3.}
  \item{c0.5A.atom }{atom indices of core positions with a total volume
    under 0.5 Angstrom\^3.}
  \item{c0.5A.xyz }{xyz indices of core positions with a total volume
    under 0.5 Angstrom\^3.}
  \item{c0.5A.resno }{residue numbers of core positions with a total volume
    under 0.5 Angstrom\^3.}
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.

  Gerstein and Altman (1995) \emph{J. Mol. Biol.}  \bold{251}, 161--175.

  Gerstein and Chothia (1991) \emph{J. Mol. Biol.} \bold{220}, 133--149.
}
\note{
  The relevance of the \sQuote{core positions} identified by this
  procedure is dependent upon the number of input structures and their
  diversity.
}
\author{ Barry Grant }
\seealso{ \code{\link{read.fasta.pdb}}, \code{\link{plot.core}},
  \code{\link{fit.xyz}} }
\examples{
\dontrun{
##--  Read kinesin alignment and respective PDB structures
aln <- read.fasta(system.file("examples/kinesin_xray.fa",package="bio3d"))
pdb.path = system.file("examples",package="bio3d")
pdbs <- read.fasta.pdb(aln, pdb.path = pdb.path, pdbext = ".ent")
}

##-- Or read previously saved kinesin data
data(kinesin)
attach(kinesin)

## Raw RMSD before superposition
gaps <- gap.inspect(pdbs$xyz)
rmsd( pdbs$xyz[,gaps$f.inds] )

## RMSD after superposition on all positions
#rmsd(pdbs$xyz[,gaps$f.inds],fit=TRUE)

## Run core.find
core <- core.find(pdbs,
                  #write.pdbs = TRUE,
                  verbose=TRUE)


## Plot volume vs length
plot(core)

## Print 0.5A^3 core and store indices
inds <- print(core, vol=0.5)

## Fit structures onto first structure based on core indices (inds$xyz)
xyz <- fit.xyz( fixed = pdbs$xyz[1,],
                mobile = pdbs,
                fixed.inds  = inds$xyz,
                mobile.inds = inds$xyz)


# RMSD after superposition on 'core' positions
rmsd( xyz[,gaps$f.inds] )


\dontrun{
# Fit structures and write out 'full' structures
xyz <- fit.xyz( fixed = pdbs$xyz[1,],
                mobile = pdbs,
                fixed.inds  = core$c0.5A.xyz,
                mobile.inds = core$c0.5A.xyz,
                pdb.path = system.file("examples/",package="bio3d"),
                pdbext = ".ent",
                outpath = "fitlsq/",
                full.pdbs = TRUE)

gaps  <- unique(which( is.na(xyz),arr.ind=TRUE )[,2])

# core fitted RMSD
rmsd(xyz[1,-gaps], xyz[,-gaps])

# original RMSD
rmsd(xyz[1,-gaps], xyz[,-gaps], fit=TRUE)


##-- Try core.find() on a trajectory
trtfile <- system.file("examples/hivp.dcd", package="bio3d")
trj <- read.dcd(trtfile)

## Read the starting PDB file to determine atom correspondence
pdbfile <- system.file("examples/hivp.pdb", package="bio3d")
pdb <- read.pdb(pdbfile)

## select calpha coords from a manageable number of frames
ca.ind <- atom.select(pdb, "calpha")$xyz
frames <- seq(1, nrow(trj), by=10)

core <- core.find( trj[frames, ca.ind], write.pdbs=TRUE )

## have a look at the various cores "vmd -m core_pruned/*.pdb"

## Lets use a 6A^3 core cutoff
inds <- print(core, vol=6)
write.pdb(xyz=pdb$xyz[inds$xyz],resno=pdb$atom[inds$atom,"resno"], file="core.pdb")


##- Fit trj onto starting structure based on core indices
xyz <- fit.xyz( fixed = pdb$xyz,
               mobile = trj,
               fixed.inds  = inds$xyz,
               mobile.inds = inds$xyz)

##write.pdb(pdb=pdb, xyz=xyz, file="new_trj.pdb")
##write.ncdf(xyz, "new_trj.nc")
}

}
\keyword{ utilities }
