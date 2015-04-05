## obtain pdb_seqres.txt from
## ftp://ftp.rcsb.org/pub/pdb/derived_data/pdb_seqres.txt
## include in data/

phmmer_local <- function(seq, exefile="phmmer", db="data/pdb_seqres.txt") {

  os1 <- .Platform$OS.type
  status <- system(exefile, ignore.stderr = TRUE, ignore.stdout = TRUE)
  
  if(!(status %in% c(0,1)))
    stop(paste("Launching external program failed\n",
               "  make sure '", exefile, "' is in your search path", sep=""))

  if(!file.exists(db))
    stop("local database not found")

  infile <- tempfile()
  outfile1 <- tempfile()
  outfile2 <- tempfile()  
  write.fasta(seq, file=infile)
  
  system(paste(exefile, "--noali --tblout", outfile2, "-o", outfile1, infile, db))
  
  return(.read_hmmerout(outfile2))
}


.read_hmmerout <- function(file) {

  trim <- function(s) {
    ##- Remove leading and trailing spaces from character strings
    s <- sub("^ +", "", s)
    s <- sub(" +$", "", s)
    s[(s=="")]<-""
    s
  }
  
  rawlines <- readLines(file)
  inds <- grep("#", rawlines)
  lines <- rawlines[-inds]
  
  acc <- trim(substr(lines, 1, 6))
  evalue <- as.numeric(trim(substr(lines, 65, 73)))
  score <- as.numeric(trim(substr(lines, 75, 80)))
  desc <- trim(substr(lines, 169, 220))
  return(data.frame(acc=acc, evalue=evalue, score=score, desc=desc,
                    stringsAsFactors=FALSE))
  
}
