seq2aln <-
function(seq2add, aln, id="seq", exefile = "muscle", file = "aln.fa") {
  ##- Add a sequence 'seq2add' to an existing alignment 'aln'
  ##  Adds at the bottom of alignment

  os1   <- .Platform$OS.type

  basealn <- tempfile()
  toaln <- tempfile()

  write.fasta(aln, file=basealn)
  if(is.vector(seq2add)) {
    write.fasta( list( id=id ,ali=seq2add), file=toaln)
  } else {
    if(is.list(seq2add)) {
      if(nrow(seq2add$ali) > 1)
        warning("Multiple sequences in 'seq2add' should be pre-aligned")
      write.fasta(seq2add, file=toaln)
    }
  }
 
  cmd <- paste(exefile, " -profile -in1 ",
               basealn, " -in2 ", toaln,
               " -out ", file, sep="")
  cat(cmd)
  
  if (os1 == "windows") {
#    system(shQuote(cmd))
    shell(shQuote(cmd))
  } else {
    system(cmd)
  }
  
  naln <- read.fasta(file)
  unlink(c(basealn, toaln))
  return(naln)
}

