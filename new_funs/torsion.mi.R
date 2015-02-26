torsion.mi <- function(xyz, pdb, normalize = TRUE, 
   xyz2tor = TRUE, ncore=NULL) {

  ncore <- setup.ncore(ncore, bigmem = TRUE)
  prev.warn <- getOption("warn")

  trj <- xyz 
  oops <- require(entropy)
  if(!oops){
    stop("Missing entropy package. Please see install.packages")
  }

  if((!is.numeric(trj)) | (!is.matrix(trj))){
    stop("'trj' object must be a numeric matrix")
  }

  if(xyz2tor) {
     cat("Calculating side-chain torsion angles...\n")
     tor <- xyz2torsion(pdb, trj, tbl = "chi1", ncore=ncore)
  } else {
     tor <- trj
  }

  cat("\nCalculating mutual information correlation...\n")

  ## Calculation of the rotamer state for each residue
  resno <- pdb$atom[pdb$calpha, "resno"]
  tor.resno <- sub("\\..*$", "", colnames(tor))
  res.inds <- which(resno %in% tor.resno)
  prot.seq <- pdbseq(pdb)
  tor.seq <- prot.seq[res.inds]

  # 0: trans; 1: gauche+; 2: gauche-
  rotamer.matrix <- matrix(0, nrow(tor), ncol(tor))
  rotamer.matrix[tor>= 0 & tor < 120] <- 1
  rotamer.matrix[tor>= -120 & tor < 0] <- 2
  # for Proline, 1: gauche+; 2: gauche-
  pro.inds <- which(tor.seq == "P")
  ii <- (tor[, pro.inds] >= 0 & tor[, pro.inds] < 90) |
        (tor[, pro.inds] >= -180 & tor[, pro.inds] < -90)
  rotamer.matrix[, pro.inds] <- 2
  rotamer.matrix[, pro.inds][ii] <- 1

  # Check state sampling and print warnings
  ns <- apply(rotamer.matrix, 2, function(x) length(unique(x)))
  if(any(ns[-pro.inds] < 3) || any(ns[pro.inds] < 2))
     warning("Some state(s) are not sampled!")

  ## Mutual information calculation
  MI.estimate <- matrix(0, nrow=length(resno), ncol=length(resno)) # Full matrix
  diag(MI.estimate) <- 1
  mi.inds <- pairwise(ncol(rotamer.matrix))
  
  if(ncore > 1) {
     mylapply <- mclapply
     options(warn = 1)
     iipb <- big.matrix(1, nrow(mi.inds), init=NA)
  } else {
     mylapply <- lapply
  } 
  pb <- txtProgressBar(min=0, max=nrow(mi.inds), style=3)
  mi <- mylapply(1:nrow(mi.inds), function(j) {
     numBins1 = numBins2 = 3
     if(mi.inds[j, 1] %in% pro.inds) numBins1 = 2
     if(mi.inds[j, 2] %in% pro.inds) numBins2 = 2
     counts <- discretize2d(
                  rotamer.matrix[, mi.inds[j,1]],
                  rotamer.matrix[, mi.inds[j,2]],
                  numBins1 = numBins1,
                  numBins2 = numBins2,
                  r1=c(0, numBins1 - 1),
                  r2=c(0, numBins2 - 1)
                  )
     suppressWarnings(mi <- mi.plugin(counts))
     if(normalize) {
        D <- (entropy::entropy(rowSums(counts)) + entropy::entropy(colSums(counts)))/2
        if(D>0) mi <- mi / D
     }

     pbj <- j
     if(ncore > 1) {
        iipb[1, j] <- 1
        pbj <- sum(!is.na(iipb[1, ]))
     } 
     setTxtProgressBar(pb, pbj)
     return(mi)
  } )
  if(ncore > 1) {
     options(warn = prev.warn)  
     gc()
  }
  close(pb)

  sub.inds <- !prot.seq %in% c("A", "G")
  sub.mi <- MI.estimate[sub.inds, sub.inds]
  sub.mi[lower.tri(sub.mi)] <- unlist(mi)
  sub.mi[upper.tri(sub.mi)] <- t(sub.mi)[upper.tri(sub.mi)]
  MI.estimate[sub.inds, sub.inds] <- sub.mi
  class(MI.estimate) <- c("dccm","matrix")
  
  return(MI.estimate)
}
