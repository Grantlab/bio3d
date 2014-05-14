"mustang" <- function(files, exefile="mustang", outfile="aln.mustang.fa",
                      verbose=TRUE) {
  ## Check if the program is executable
  os1 <- .Platform$OS.type
  status <- system(paste(exefile, "--version"),
                   ignore.stderr = TRUE, ignore.stdout = TRUE)
  
  if(!(status %in% c(0,1)))
    stop(paste("Launching external program failed\n",
               "  make sure '", exefile, "' is in your search path", sep=""))
  
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
    cat("Running command\n", cmd, "\n")
  
  if (os1 == "windows")
    success <- shell(shQuote(cmd), ignore.stderr = !verbose, ignore.stdout = !verbose)
  else
    success <- system(cmd, ignore.stderr = !verbose, ignore.stdout = !verbose)

  if(success!=0)
    stop(paste("An error occurred while running command\n '",
               exefile, "'", sep=""))
  
  aln <- read.fasta(paste(outfile, ".afasta", sep=""))
  rownames(aln$ali) <- paste(dirn, rownames(aln$ali), sep="/")
  aln$id <- rownames(aln$ali)

  unlink(infile); unlink(outfile);

  if(!is.null(outfile))
    write.fasta(aln, file=outfile)
  
  return(aln)
}
