###
### Examples from NMA Vignette
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

#######################################
## Basic usage                        #
#######################################

### Read PDB and Calculate Normal Modes
pdb <- read.pdb("1hel")
modes <- nma(pdb)

pause()

### Print a summary
print(modes)

pause()

### Plot the nma object for a quick overview
plot(modes)

pause()

### Calculate cross-correlations
cm <- dccm(modes)

pause()

### Plot correlation map
plot(cm, sse=pdb.open)

pause()

### Calculate modes with force field ANM
modes.anm <- nma(pdb, ff="anm")

pause()

### Investigate modes similarity with RMSIP
r <- rmsip(modes, modes.anm)

pause()

### Plot RMSIP results
plot(r, xlab="ANM", ylab="C-alpha FF")

pause()


#######################################
## Ensemble NMA                       #
#######################################


### Download a set of Kinase structures
ids <- c("4b7t_A", "2exm_A", "1opj_A", 
         "4jaj_A", "1a9u_A", "1tki_A")

### Download and split by chain ID
raw.files <- get.pdb(ids, path="raw_pdbs")

pause()

### Split PDB files by chain ID
files     <- pdbsplit( raw.files, ids )

pause()

### Align structures
pdbs <- pdbaln(files)

pause()

### View sequence identity
summary( c(seqidentity(pdbs)) )

pause()

### Calculate modes of aligned proteins
modes <- nma.pdbs(pdbs, fit=TRUE)

pause()

### Plot fluctuations
plot(modes, pdbs, type="h")

pause()
