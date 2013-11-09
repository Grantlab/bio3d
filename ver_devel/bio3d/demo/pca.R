###
### Example of PCA on a collection of PKA structures
### Authors Lars Skjaerven
###         Xin-Qiu Yao
###         Barry J Grant
###
require(bio3d); require(graphics);

pause <- function() {
  cat("Press ENTER/RETURN/NEWLINE to continue.")
  readLines(n=1)
  invisible()
}

#############################################
## Basic PCA of related X-ray structures    #
#############################################

## Specify PDB identifiers
ids <- c("1cdk_A", "3agm_A", "1cmk_E",
         "3dnd_A", "1q8w_A")

## Download PDBs
raw.files <- get.pdb(ids)

## Split PDBs by chain ID
files <- pdbsplit(raw.files, ids)

pause()

## Sequence/structure alignment
pdbs <- pdbaln(files)

pause()

## Find invariant core
core <- core.find(pdbs)

pause()

## Fit structures to core region
xyz <- fit.xyz(fixed=pdbs$xyz[1,], mobile=pdbs,
               fixed.inds=core$c1A.xyz, mobile.inds=core$c1A.xyz)
               ##outpath="core_fit/", full.pdbs=T, het2atom=T)

pause()

## Inspect gaps
gaps.pos <- gap.inspect(pdbs$xyz)

## Perform PCA on non-gap containing positions
pc.xray <- pca.xyz(xyz[,gaps.pos$f.inds])

pause()

## Plot x-ray results
plot(pc.xray)

