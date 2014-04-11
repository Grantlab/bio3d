"mustang" <- function(files, exefile="mustang", verbose=TRUE) {
  if(!all(file.exists(files)))
    stop(paste("Missing files:", paste(files[ !file.exists(files) ], collapse=", ")))

  infile <- tempfile()
  outfile <- tempfile()
  dirn <- unique(dirname(files))

  if(length(dirn)>1)
    stop("All files must be in one directory")
  
  files <- basename(files)

  rawlines <- NULL
  rawlines <- c(rawlines, paste(">", dirn))

  for ( i in 1:length(files) )
    rawlines <- c(rawlines, paste("+", files[i], sep=""))
  
  write.table(rawlines, file=infile, quote=FALSE,
              row.names=FALSE, col.names=FALSE)
  
  cmd <- paste(exefile, "-f", infile, "-o", outfile, "-F fasta")

  if(verbose)
    print(cmd)
  system(cmd)
  
  aln <- read.fasta(paste(outfile, ".afasta", sep=""))
  rownames(aln$ali) <- paste(dirn, rownames(aln$ali), sep="/")
  aln$id <- rownames(aln$ali)

  ## check for only NA columns
  #gaps <- gap.inspect(aln$ali)
  #inds <- which(gaps$col<ncol(aln$ali))
  #aln$ali=aln$ali[,inds]

  ##pdbs <- read.fasta.pdb(aln)
  return(aln)
}
