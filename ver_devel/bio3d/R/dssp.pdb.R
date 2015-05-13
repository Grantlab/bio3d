## NOTE:
##   We do not support old-version DSSP any longer
##   Please update your DSSP program to the newest version
"dssp.pdb" <-
  function (pdb, exefile = "dssp", resno=TRUE, full=FALSE, verbose=FALSE, ...) {

    ## Log the call
    cl <- match.call()
    
    ## Check if the program is executable
    os1 <- .Platform$OS.type
    status <- system(paste(exefile, "--version"),
                     ignore.stderr = TRUE, ignore.stdout = TRUE)
    
###    if(!(status %in% c(0,1)))
###      stop(paste("Launching external program 'DSSP' failed\n",
###                 "  make sure '", exefile, "' is in your search path", sep=""))
    
    ## check atom composition - need backbone atoms to continue SSE analysis
    checkatoms <- TRUE
    if(checkatoms) {
      inds <- atom.select(pdb, "backbone", verbose=verbose)
      tmp <- trim.pdb(pdb, inds)
      
      resid <- paste(tmp$atom$resno, tmp$atom$chain, sep="-")
      musthave <- c("C", "CA", "N", "O")

      incomplete <- sapply(unique(resid), function(x) {
        inds <- which(resid==x)
        elety <- sort(tmp$atom$elety[inds])
        if(!all(musthave %in% elety))
          return(TRUE)
        else
          return(FALSE)
      })

      if(all(incomplete))
        stop("No residues found with a complete set of backbone atoms")
      if(any(incomplete))
        warning(paste("Residues with missing backbone atoms detected:",
                      paste(unique(resid)[incomplete], collapse=", "), 
                      collapse=" "))
    }
    
    infile <- tempfile()
    outfile <- tempfile()
    write.pdb(pdb, file = infile)
    cmd <- paste(exefile, infile, outfile)

    if(verbose)
      cat(paste("Running command:\n ", cmd , "\n"))
    
    if(os1 == "windows")
      success <- shell(cmd, 
                       ignore.stderr = !verbose, ignore.stdout = !verbose)
    else 
      success <- system(cmd,
                        ignore.stderr = !verbose, ignore.stdout = !verbose)
    
    if(success!=0)
      stop(paste("An error occurred while running command\n '",
                 cmd, "'", sep=""))
   
    out <- read.dssp(outfile, resno=resno, full=full)
    out$call <- cl
    unlink(c(infile, outfile))

    return(out)
}
