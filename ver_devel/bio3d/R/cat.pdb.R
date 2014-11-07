cat.pdb <- function(pdb1, pdb2, renumber=TRUE, rechain=TRUE) {
  if(missing(pdb1) | !is.pdb(pdb1) | missing(pdb2) | !is.pdb(pdb2) )
    stop("provide 'pdb1' and 'pdb2' as obtained from read.pdb()")
  
  cl <- match.call()
  
  if(rechain) {
    j <- 1
    chains <- unique(pdb1$atom$chain)
    for(i in 1:length(chains)) {
      inds <- which(pdb1$atom$chain==chains[i])
      pdb1$atom$chain[inds] <- LETTERS[j]
      j <- j+1
    }
    
    chains <- unique(pdb2$atom$chain)
    for(i in 1:length(chains)) {
      inds <- which(pdb2$atom$chain==chains[i])
      pdb2$atom$chain[inds] <- LETTERS[j]
      j <- j+1
    }
  }
  
  new <- pdb1
  new$atom       <- rbind(pdb1$atom, pdb2$atom)
  new$xyz        <- cbind(pdb1$xyz, pdb2$xyz)
  new$atom$eleno <- seq(1, nrow(new$atom))
  
  if(renumber)
    prev.resno <- 0
  else
    prev.resno <- new$atom$resno[1]

  prev.resid   <- new$atom$resid[1]
  prev.chain   <- new$atom$chain[1]
  new.resno    <- prev.resno
  
  for(i in 1:nrow(new$atom)) {
    now.resno  <- new$atom$resno[i]
    now.resid  <- new$atom$resid[i]
    now.chain  <- new$atom$chain[i]

    if(rechain & (now.chain!=prev.chain))
      new.resno <- 0
    
    if(( now.resno!=prev.resno) | (now.resid!=prev.resid)) {
      new.resno <- new.resno+1
    }
    new$atom$resno[i] <- new.resno

    prev.resno <- now.resno
    prev.resid <- now.resid
    prev.chain <- now.chain
  }
  
  new$sheet <- NULL
  new$helix <- NULL
  new$seqres <- NULL
  ca.inds <- atom.select(new, "calpha", verbose=FALSE)
  new$calpha <- seq(1, nrow(new$atom)) %in% ca.inds$atom
  new$call <- cl
  

  ## check connectivity
  chains <- unique(new$atom$chain)
  for(i in 1:length(chains)) {
    sele <- atom.select(new, chain=chains[i], verbose=FALSE)
    tmp <- trim.pdb(new, sele)
    
    if(!inspect.connectivity(tmp))
      warning(paste("possible chain break in molecule: chain", chains[i]))
  }
  
  return(new)
}
