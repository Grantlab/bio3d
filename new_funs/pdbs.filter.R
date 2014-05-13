## Use for trimming a 3dalign object, either by removing structures,
## or by removing columns 

pdbs.filter <- function(pdbs, row.inds=NULL, col.inds=NULL) {
  if(!inherits(pdbs, "3dalign"))
    stop("input 'pdbs' should be a list object as obtained from 'read.fasta.pdb'")
  
  ## Log the call
  cl <- match.call()
  
  if(is.null(col.inds))
    col.inds <- seq(1, ncol(pdbs$resno), by=1)
  if(is.null(row.inds))
    row.inds <- seq(1, nrow(pdbs$resno), by=1)
  
  gaps <- gap.inspect(pdbs$resno[row.inds, col.inds])
  col.inds <- which(gaps$col < length(row.inds))
  
  new <- NULL
  new$id    =pdbs$id[row.inds]
  new$xyz   =pdbs$xyz[row.inds, atom2xyz(col.inds)]
  new$resno =pdbs$resno[row.inds, col.inds]
  new$b     =pdbs$b[row.inds, col.inds]
  new$chain =pdbs$chain[row.inds, col.inds]
  new$ali   =pdbs$ali[row.inds, col.inds]
  new$resid =pdbs$resid[row.inds, col.inds]
  new$call  =cl
  
  class(new) <- c("3dalign", "fasta")
  return(new)
}
