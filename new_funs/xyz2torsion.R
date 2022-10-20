"xyz2torsion" <-
function(pdb, xyz, tbl = c("basic", "mainchain", "sidechain", 
    "all", "phi", "psi", paste("chi", 1:5, sep="")), ncore = NULL) {

  if(ncol(pdb$xyz) != ncol(xyz)) 
     stop("Number of atoms in PDB doesn't match xyz")
  ncore <- setup.ncore(ncore, bigmem = TRUE)
  
  tor.names <- c("phi", "psi", paste("chi", 1:5, sep=""))
  tbl <- match.arg(tbl)
  tbl <- switch(tbl, 
     "basic" = tor.names[1:3],
     "mainchain" = tor.names[1:2], 
     "sidechain" = tor.names[3:7],
     "all" = tor.names,
     tbl
  )
  tor <- torsion.pdb(pdb)
  fill <- !is.na(t(tor$tbl[, tbl]))
  resno <- rep(pdb$atom[pdb$calpha, "resno"], each = length(tbl))
  cname <- paste(resno[fill], rep(tbl, ncol(fill))[fill], sep=".")
  n <- sum(fill)

  pb <- txtProgressBar(min=0, max=nrow(xyz), style=3)
  if(ncore > 1) {
     tor <- big.matrix(nrow(xyz), n, init=NA, type="double")
     iipb <- big.matrix(1, nrow(xyz), init=NA)
     mclapply(1:nrow(xyz), function(i) {
        pdb$xyz <- xyz[i, ]
        tor[i, ] <- t(torsion.pdb(pdb)$tbl[, tbl])[fill]
        iipb[1, i] <- 1
        setTxtProgressBar(pb, sum(!is.na(iipb[1,])))
        return()
     } )
     tor <- tor[,]; gc()
  } else {
     tor <- t( sapply(1:nrow(xyz), function(i) {
        pdb$xyz <- xyz[i, ]
        tor <- t(torsion.pdb(pdb)$tbl[, tbl])[fill]
        setTxtProgressBar(pb, i)
        return(tor)
     } ) )
  }
  close(pb)
  colnames(tor) <- cname
  return(tor)
}
