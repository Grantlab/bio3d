cat.pdb <- function(..., renumber=TRUE, rechain=TRUE) {
  cl <- match.call()
    
  objs <- list(...)
  are.null <- unlist(lapply(objs, is.null))
  objs <- objs[!are.null]

  if(length(objs)<2)
    stop("provide multiple (>1) PDB objects")

  if(any(!unlist(lapply(objs, is.pdb))))
    stop("provide PDB objects as obtained from read.pdb()")
  
  rechainfun <- function(x, j) {
    chains <- unique(x$atom$chain)
    for(i in 1:length(chains)) {
      inds <- which(x$atom$chain==chains[i])
      x$atom$chain[inds] <- LETTERS[j]
      j <- j+1
    }
    return(list(pdb=x, j=j))
  }

  if(rechain) {
    j <- 1
    for(i in 1:length(objs)) {
      tmp <- rechainfun(objs[[i]], i)
      objs[[i]] <- tmp$pdb
      j <- tmp$j
    }
  }
  
  new <- objs[[1]]
  for(i in 2:length(objs)) {
    new$atom     <- rbind(new$atom, objs[[i]]$atom)
    new$xyz      <- cbind(new$xyz, objs[[i]]$xyz)
  }
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
