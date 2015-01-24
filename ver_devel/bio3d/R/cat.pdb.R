cat.pdb <- function(..., renumber=FALSE, rechain=FALSE) {
  cl <- match.call()
    
  objs <- list(...)
  are.null <- unlist(lapply(objs, is.null))
  objs <- objs[!are.null]

  if(length(objs)<1) return(NULL)

  if(any(!unlist(lapply(objs, is.pdb))))
    stop("provide PDB objects as obtained from read.pdb()")
  
  ## avoid NA as chain ID
  na.inds <- lapply(objs, function(x) is.na(x$atom$chain))
  if(any(unlist(na.inds))) {
    na.inds <- which(unlist(lapply(na.inds, function(x) any(x))))
    for(i in na.inds) {
      tmp <- objs[[i]]
      tmp$atom$chain[ is.na(tmp$atom$chain) ] <- " "
      objs[[i]] <- tmp
    }
  }

  ## assign new chain identifiers
  if(rechain) {
    k <- 1
    for(i in 1:length(objs)) {
      x <- objs[[i]]
      chains <- unique(x$atom$chain)
      for(j in 1:length(chains)) {
        inds <- which(objs[[i]]$atom$chain==chains[j])
        x$atom$chain[inds] <- LETTERS[k]
        if(!is.null(x$helix)) x$helix$chain[] <- LETTERS[k]
        if(!is.null(x$sheet)) x$sheet$chain[] <- LETTERS[k]
        k <- k+1
      }
      objs[[i]] <- x
    }
  }

  ## concat objects
  new <- objs[[1]]
  if(length(objs) > 1) { 
     for(i in 2:length(objs)) {
       new$atom     <- rbind(new$atom, objs[[i]]$atom)
       new$xyz      <- cbind(new$xyz, objs[[i]]$xyz)
       new$seqres <- c(new$seqres, objs[[i]]$seqres)
       new$helix <- c(new$helix, objs[[i]]$helix)
       new$sheet <- c(new$sheet, objs[[i]]$sheet)
     }
  }
  ## merge SSE info
  for(i in c("helix", "sheet")) {
     sse <- new[[i]]
     if(!is.null(sse)) {
        coms <- names(sse)
        names(sse) <- NULL # avoid nested naming in results
        inds <- which(!duplicated(coms))
        for(j in inds) 
           sse[[j]] <- do.call(c, sse[coms %in% coms[j]])
        sse <- sse[inds]
        names(sse) <- coms[inds]
        new[[i]] <- sse
     }
  }

  ## renumber residues
  new <- try(clean.pdb(new, consecutive = !rechain, 
                force.renumber = renumber, verbose=FALSE))
  if(inherits(new, "try-error")) 
     stop("cat.pdb(): Bad format pdb generated. Try rechain=TRUE and/or renumber=TRUE")

  ## build new PDB object
  new$call <- cl

  ## remap " " chain IDs to NA values
  new$atom$chain[ new$atom$chain==" " ] <- as.character(NA)

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
