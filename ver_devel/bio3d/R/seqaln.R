"seqaln" <-
function(aln, id=NULL,
                   exefile = "muscle",
                   outfile = "aln.fa",
                   protein = TRUE,
                   seqgroup = FALSE,
                   refine = FALSE,
                   extra.args = "",
                   verbose = FALSE) {

  ## Check if the program is executable
  os1 <- .Platform$OS.type
  status <- system(paste(exefile, "-version"),
                   ignore.stderr = TRUE, ignore.stdout = TRUE)
  
  if(!(status %in% c(0,1)))
    stop(paste("Launching external program 'MUSCLE' failed\n",
               "  make sure '", exefile, "' is in your search path", sep=""))
  
  as.aln <- function(mat, id=NULL) {
    if(is.null(id))
      id=paste("seq",1:nrow(mat),sep="")
    return(list(id=id, ali=mat))
  }

  if( (!is.list(aln)) | is.na(aln['id']) )
    aln<-as.aln(aln,id=id)


  toaln <- tempfile()  
  write.fasta(aln, file=toaln)
  
  if(is.null(outfile)) fa <- tempfile() 
  else fa <- outfile

###  if(!seqgroup)  extra.args <- paste(extra.args,"-stable")
  if(refine) extra.args <- paste(extra.args,"-refine")
  if(protein) {
    extra.args <- paste(extra.args,"-seqtype protein")
  } else { extra.args <- paste(extra.args,"-seqtype dna") }
   
  cmd <- paste(exefile, " -in ",toaln," -out ",
               fa," ",extra.args, sep="")
  if(verbose)
    cat(paste("Running command:\n ", cmd , "\n"))
  
  ## Run command
  if (os1 == "windows")
    success <- shell(shQuote(cmd), ignore.stderr = !verbose, ignore.stdout = !verbose)
  else
    success <- system(cmd, ignore.stderr = !verbose, ignore.stdout = !verbose)
  
  if(success!=0)
    stop(paste("An error occurred while running command\n '",
               exefile, "'", sep=""))

  ### Update for muscle v3.8 with no "-stable" option
  ###  Thu Aug 26 18:29:38 PDT 2010
  naln <- read.fasta(fa, rm.dup=FALSE)
  if(!seqgroup) {
    ord <- match(aln$id, naln$id)
##    if( any( duplicated(ord) ) ) {
##      stop(" Duplicated sequence id's, so can't perserve input order\n\tPlease run with 'seqgroup=TRUE' option")
##    }
    naln$id <- naln$id[ord]
    naln$ali <- naln$ali[ord,]
  }
  ####
  
  unlink(toaln)
  if(is.null(outfile)) unlink(fa)
  return(naln)
}

