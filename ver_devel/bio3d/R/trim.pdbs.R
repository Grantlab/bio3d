## Use for trimming a pdbs object, either by removing structures,
## or by removing columns 

trim.pdbs <- function(pdbs, row.inds=NULL, col.inds=NULL, ...) {
  if(!inherits(pdbs, "pdbs"))
    stop("input 'pdbs' should be a list object as obtained from 'read.fasta.pdb'")
  
  ## Log the call
  cl <- match.call()
  
  if(is.null(row.inds))
    row.inds <- seq(1, nrow(pdbs$resno), by=1)
  if(is.null(col.inds)) {
    gaps <- gap.inspect(pdbs$resno[row.inds,,drop=FALSE])
    col.inds <- which(gaps$col < dim(pdbs$resno[row.inds,,drop=FALSE])[1L])
  }
  
  if(any(col.inds<0))
    col.inds.xyz <- atom2xyz(abs(col.inds)) * sign(rep(col.inds, each=3))
  else
    col.inds.xyz <- atom2xyz(col.inds)
  
  new <- NULL
  new$id    =pdbs$id[row.inds]
  new$xyz   = as.xyz(pdbs$xyz[row.inds, col.inds.xyz, drop=FALSE])
  new$resno =pdbs$resno[row.inds, col.inds, drop=FALSE]
  new$b     =pdbs$b[row.inds, col.inds, drop=FALSE]
  new$chain =pdbs$chain[row.inds, col.inds, drop=FALSE]
  new$ali   =pdbs$ali[row.inds, col.inds, drop=FALSE]
  new$resid =pdbs$resid[row.inds, col.inds, drop=FALSE]
  new$sse   =pdbs$sse[row.inds, col.inds, drop=FALSE]
  new$call  =cl
 
  if(!is.null(pdbs$all)) {
    col.inds.all <- which(pdbs$all.grpby %in% abs(col.inds))
    col.inds.all <- col.inds.all * sign(rep(col.inds, rle(pdbs$all.grpby[col.inds.all])$length))

    if(any(col.inds.all<0))
      col.inds.all.xyz <- atom2xyz(abs(col.inds.all)) * sign(rep(col.inds.all, each=3))
    else
      col.inds.all.xyz <- atom2xyz(col.inds.all)

    new$all       = as.xyz(pdbs$all[row.inds, col.inds.all.xyz, drop=FALSE])
    new$all.elety = pdbs$all.elety[row.inds, col.inds.all, drop=FALSE]
    new$all.resno = pdbs$all.resno[row.inds, col.inds.all, drop=FALSE]
    new$all.resid = pdbs$all.resid[row.inds, col.inds.all, drop=FALSE]
    new$all.grpby   = vec2resno(1:ncol(new$ali), pdbs$all.grpby[col.inds.all])
  }

  class(new) <- c("pdbs", "fasta")
  return(new)
}
