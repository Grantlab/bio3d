"dssp.pdbs" <- function(pdbs) {

  gaps.res <- gap.inspect(pdbs$ali)
  ##gaps.pos <- gap.inspect(pdbs$xyz)

  sse <- matrix(NA, ncol=ncol(pdbs$resno), nrow=nrow(pdbs$resno))
  
  for ( i in 1:length(pdbs$id) ) {
    if(!file.exists(pdbs$id[i]))
      stop(paste(pdbs$id[i], "does not exist"))
    
    tmp.pdb = read.pdb(pdbs$id[i])
    tmp.sse = dssp(tmp.pdb, resno=FALSE, full=FALSE, verbose=FALSE)
    
    tmp.sse$sse[ tmp.sse$sse == " " ] = "-"
    sse[i,  which(gaps.res$bin[i,]==0)] = tmp.sse$sse
  }

  sse[is.na(sse)]="-"
  return(sse)
}
