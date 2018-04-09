`chain.pdb` <-
function(pdb, ca.dist=4, bond=TRUE, bond.dist=1.5, blank="X", rtn.vec=TRUE) {
  ##- Find possible chian breaks
  ##  i.e. Concetive Caplpa's that are further than 'ca.dist' apart,
  ##        print basic chain info and rtn a vector of chain ids
  ##        consisting of the 26 upper-case letters of the Roman
  ##        alphabet
  ##
  ## chn <- chain.pdb(pdb)
  ## pdb$atom[,"chain"] <- chain.pdb(pdb)
  ##

  ## Distance between concetive C-alphas
  ca <- atom.select(pdb, "calpha", verbose=FALSE)
  if(length(ca$atom) <=1) {
    d <- 0
  } else {
    if(bond) {
      nn <- atom.select(pdb, "protein", elety="N", verbose=FALSE)
      cc <- atom.select(pdb, "protein", elety="C", verbose=FALSE)
      if(length(ca$atom) != length(nn$atom) ||
         length(ca$atom) != length(cc$atom)) {
        stop("Peptide bond atoms (N/C) and C-alpha atoms mismatch (try bond=FALSE).")
      } else {
        xyz1 <- pdb$xyz[cc$xyz]
        xyz2 <- pdb$xyz[nn$xyz][-c(1:3)]
        d <- dist.xyz(xyz1, xyz2, all.pairs=FALSE)
        d <- d[!is.na(d)]
      }
    } else { 
      xyz <- matrix(pdb$xyz[ca$xyz], nrow=3)
      if(length(ca$atom) == 2) 
        d <- sqrt( sum( apply(xyz , 1, diff)^2 ) )
      else
        d <- sqrt( rowSums( apply(xyz , 1, diff)^2 ) )
    }
  }

  ## Chain break distance check
  if(bond) {
    ind <- which(d > bond.dist)
  } else {
    ind <- which(d > ca.dist)
  }
  len <- diff( c(1,ind,length(d)) )

  cat(paste("\t Found",length(ind), "possible chain breaks\n"))
  if(length(ind) > 0) {
     cat(paste("\t  After resno(s):",
           paste( pdb$atom[ca$atom,"resno"][(ind)], collapse=", " ),"\n" ))
     cat(paste("\t  Chain length(s):",
           paste(len+1, collapse=", " ),"\n" ))
  }

  ## Make a chain id vector
  if(rtn.vec) {
    resno.ind <- as.numeric(c(1, sort(as.numeric(c(ind,(ind+1)))), (length(d)+1)
))
    ## Renumber residues first, in case that original resnos are not
    ## consecutive crossing multiple chains
    res <- paste(pdb$atom[, "chain"], pdb$atom[, "resno"], pdb$atom[, "insert"], sep="_")
    pdb$atom[, "resno"] <- vec2resno(1:length(unique(res)), res) 

    resno.val <- pdb$atom[ca$atom,"resno"][resno.ind]
    resno.val <- matrix(as.numeric(resno.val),nrow=2)

    vec <- rep(blank, nrow(pdb$atom))
    if(is.na(resno.val[1])) {
      return(vec)
    } else if(is.na(resno.val[2])) {
      resno.val[2] <- resno.val[1]
    }
    
    for(i in 1:(length(resno.val)/2)) {
      sel.ind <- atom.select(pdb,
                             resno=c(resno.val[1,i]:resno.val[2,i]),
                             verbose=FALSE)
      vec[sel.ind$atom]=LETTERS[i]
    }
    return(vec)
  }
}

