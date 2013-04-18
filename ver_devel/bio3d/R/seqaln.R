"seqaln" <-
function(aln, id=NULL,
                   exepath = "",
                   file = "aln.fa",
                   protein = TRUE,
                   seqgroup = FALSE,
                   refine = FALSE,
                   extra.args = "") {


  as.aln <- function(mat, id=NULL) {
    if(is.null(id))
      id=paste("seq",1:nrow(mat),sep="")
    return(list(id=id, ali=mat))
  }

  if( (!is.list(aln)) | is.na(aln['id']) )
    aln<-as.aln(aln,id=id)


  os1   <- .Platform$OS.type
  toaln <- tempfile()  
  write.fasta(aln, file=toaln)
  
###  if(!seqgroup)  extra.args <- paste(extra.args,"-stable")
  if(refine) extra.args <- paste(extra.args,"-refine")
  if(protein) {
    extra.args <- paste(extra.args,"-seqtype protein")
  } else { extra.args <- paste(extra.args,"-seqtype dna") }
  
  cmd <- paste(exepath, "muscle -in ",toaln," -out ",
               file," ",extra.args, sep="")
  cat(cmd)
  
  if (os1 == "windows") {
#    system(shQuote(cmd))
    shell(shQuote(cmd))
  } else {
    system(cmd)
  }

  ### Update for muscle v3.8 with no "-stable" option
  ###  Thu Aug 26 18:29:38 PDT 2010
  naln <- read.fasta(file, rm.dup=FALSE)
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
  return(naln)
}

