phmmer_local <- function(seq) {
  exefile <- configuration$hmmer$exefile
  db <- configuration$hmmer$pdbseq
  
  status <- system(exefile, ignore.stderr = TRUE, ignore.stdout = TRUE)
  if(!(status %in% c(0,1))) {
    stop(paste("Launching external program failed\n",
               "  make sure '", exefile, "' is in your search path", sep=""))
  }
  
  if(!file.exists(db)) {
    stop("local database not found")
  }

  infile <- tempfile()
  outfile1 <- tempfile()
  outfile2 <- tempfile()  
  write.fasta(seq, file=infile)

  system(paste(exefile, "--noali --tblout", outfile2, "-o", outfile1, infile, db))
  
  return(.read_hmmerout(outfile2))
}


.read_hmmerout <- function(file) {

  rawlines <- readLines(file)
  inds <- grep("#", rawlines)
  lines <- rawlines[-inds]

  ## FORMAT: 4v8r_AB              -          seq1                 -            2.4e-15   60.6   5.8   3.5e-10   43.6   1.2   2.1   2   0   0   2   2   2   2 mol:protein length:527  T-COMPLEX PROTEIN 1 SUBUNIT BETA

  acc <- trim(substr(lines, 1, 7)) 
  evalue <- as.numeric(trim(substr(lines, 65, 73)))
  score <- as.numeric(trim(substr(lines, 75, 80)))
  desc <- trim(substr(lines, 169, 220))
  return(data.frame(acc=acc, evalue=evalue, score=score, desc=desc,
                    stringsAsFactors=FALSE))
  
}


