"dssp.pdbs" <- function(pdbs, ...) {
  if(!is.pdbs(pdbs))
    stop("provide a pdbs object as obtained from pdbaln()")
  
  dots <- list(...)
  if(any(c("resno", "full") %in% names(dots)))
    stop("arguments resno and full not allowed in dssp.pdbs()")
  
  gaps.res <- gap.inspect(pdbs$ali)
  sse <- matrix(NA, ncol=ncol(pdbs$resno), nrow=nrow(pdbs$resno))
  
  for ( i in 1:length(pdbs$id) ) {
    if(!file.exists(pdbs$id[i]))
      stop(paste(pdbs$id[i], "does not exist"))
    
    tmp.pdb = read.pdb(pdbs$id[i])
    tmp.sse = dssp.pdb(tmp.pdb, resno=FALSE, full=FALSE, ...)
    
    tmp.sse$sse[ tmp.sse$sse == " " ] = "-"
    sse[i,  which(gaps.res$bin[i,]==0)] = tmp.sse$sse
  }

  sse[is.na(sse)]="-"
  return(sse)
}
