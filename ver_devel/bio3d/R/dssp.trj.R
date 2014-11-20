dssp.trj <- function(pdb, trj, skip=1000, threshold=3, file.head="", ...) {
  if(!is.pdb(pdb))
    stop("provide a pdb object as obtained from function 'read.pdb'")
  
  ##Filtering
  filter.raw <- seq(1,dim(trj)[1],skip)
  
  ##DSSP calculation
  dssp.ref  <- dssp(pdb, ...)
  helix.ref <- sum(dssp.ref$helix$length)/sum(pdb$calpha)*100
  sheet.ref <- sum(dssp.ref$sheet$length)/sum(pdb$calpha)*100
  
  unfold.frames <- NULL
  for (frame in filter.raw){
    pdb.temp     <- pdb
    pdb.temp$xyz <- trj[frame,]
    dssp.test    <- dssp(pdb.temp, ...)
    helix.test   <- sum(dssp.test$helix$length)/sum(pdb$calpha)*100
    sheet.test   <- sum(dssp.test$sheet$length)/sum(pdb$calpha)*100
    helix.diff   <- helix.ref - helix.test
    sheet.diff   <- sheet.ref - sheet.test
    
    if (helix.diff >= helix.ref/threshold || sheet.diff >= sheet.ref/threshold) {
      fname <- paste(file.head, frame, '.pdb', sep="")
      cat('Possible unfolding event(s) in frame', frame, '\n')
      cat('  ... saving frame to file', fname, '\n')
      
      write.pdb(pdb=pdb, file=fname, xyz=trj[frame,])
      unfold.frames <- c(unfold.frames, frame)
    } 
  }
  
  if(is.null(unfold.frames)) {
    cat('DSSP completed with no detected unfolding frames.')
  }
  else {
    cat(length(unfold.frames), 'unfolding frame(s) detected.')
  }
  cat("\n")
  
  return(unfold.frames)
}
