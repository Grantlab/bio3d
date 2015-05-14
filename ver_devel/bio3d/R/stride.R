"stride" <-
function(pdb, exefile = "stride", resno=TRUE) {
  
  ## Log the call
  cl <- match.call()
  
  infile  <- tempfile()
  outfile <- tempfile()
  write.pdb(pdb, file=infile)

  os1 <- .Platform$OS.type
  if(os1 == "windows") {
     shell( paste(exefile," -f",outfile," ",infile,sep="") )
  } else {
     system( paste(exefile," -f",outfile," ",infile,sep="") )
  }

  out <- read.stride(outfile, resno=resno)
  out$call <- cl
  unlink(c(infile, outfile))
  
  return(out)
}

