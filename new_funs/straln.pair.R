# Pairwise structure alignment based on dynamic programming
# Can do both single align or align to align
# Assume pre-fitted: fit=TRUE not supported yet...
straln.pair <- function(pdbs1, pdbs2, gap.penalty = -9,
                          ret.score = FALSE, fit = FALSE, 
                          model = c("ca", "backbone"),
                          method = c("average", "maximal"), 
                          ncore = NULL) {
   ncore <- setup.ncore(ncore)
   model <- match.arg(model)
   method <- match.arg(method)
   
   mylapply <- lapply
   if(ncore > 1) mylapply <- mclapply

   if(is.pdb(pdbs1)) {
      pdbs1 <- list(id=pdbs1$call$file, ali=matrix(pdbseq(pdbs1), nrow=1),
              xyz=matrix(pdbs1$xyz[pdbs1$calpha], nrow=1))
   }
   if(is.pdb(pdbs2)) {
      pdbs2 <- list(id=pdbs2$call$file, ali=matrix(pdbseq(pdbs2), nrow=1),
              xyz=matrix(pdbs2$xyz[pdbs2$calpha], nrow=1))
   }
  
#   gap.penalty <- - gap.penalty^2
 
   # build scoring matrix
   nseq1 <- nrow(pdbs1$ali)
   nseq2 <- nrow(pdbs2$ali)
   np <- nseq1 * nseq2
   score.all <- mylapply(1:np, function(x) {
      i <- floor((x-1) / nseq2) + 1
      j <- (x-1) %% nseq2 + 1
#      sc <- - dist.xyz(pdbs1$xyz[i, ], pdbs2$xyz[j, ])
      sc <- switch(model,
         ca = - (dist.xyz(pdbs1$xyz[i, ], pdbs2$xyz[j, ]))^2,
         backbone = {
            xyz1 <- matrix(pdbs1$xyz[i, ], ncol = 4*3, byrow=TRUE)
            xyz2 <- matrix(pdbs2$xyz[j, ], ncol = 4*3, byrow=TRUE)
            - (dist.xyz(xyz1, xyz2))^2 / 4
         } )
      g1 <- which(is.gap(pdbs1$ali[i, ]))
      g2 <- which(is.gap(pdbs2$ali[j, ]))
      sc[g1, g2] <- 0   # no penalty on gap-gap
      sc[is.na(sc)] <- gap.penalty 
      return(as.vector(sc))
   })
   score.all <- do.call(rbind, score.all)

   # collapse 
   score <- switch(method,
       maximal = apply(score.all, 2, max, na.rm=TRUE),
       average = apply(score.all, 2, mean, na.rm=TRUE) )
   
   score.mat <- matrix(score, ncol(pdbs1$ali), ncol(pdbs2$ali))
   
   # build local optimization matrix 
   mat <- matrix(NA, nrow(score.mat)+1, ncol(score.mat)+1)
   mat[1, ] <- c(0, cumsum(rep(gap.penalty, ncol(mat)-1)))
   mat[, 1] <- c(0, cumsum(rep(gap.penalty, nrow(mat)-1)))
   arg <- matrix(NA, nrow(mat), ncol(mat))
   arg[1, ] <- c("m", rep("i", ncol(arg)-1))
   arg[, 1] <- c("m", rep("d", nrow(arg)-1))

   st <- c("m", "i", "d")
   for(j in 2:nrow(mat)) {   # ref
      for(k in 2:ncol(mat)) {  # target
         sc <- mat[j-1, k-1] + score.mat[j-1, k-1]   # match
         sc <- c(sc, mat[j, k-1] + gap.penalty)      # insert
         sc <- c(sc, mat[j-1, k] + gap.penalty)      # delete
         arg[j, k] <- st[which.max(sc)]
         mat[j, k] <- max(sc)
      }
   }

   # backward trace optimal path
   j = nrow(score.mat); k = ncol(score.mat); n=j+k
   ali <- matrix("-", nseq1+nseq2, n)
   xyz <- matrix(NA,  nseq1+nseq2, n*3)
   while(j > 0 || k > 0) {
      switch(arg[j+1, k+1],
         m = {
            ali[, n] <- c(pdbs1$ali[, j], pdbs2$ali[, k])
            xyz[, atom2xyz(n)] <- rbind(pdbs1$xyz[, atom2xyz(j)], pdbs2$xyz[, atom2xyz(k)])
            j = j-1; k = k-1
         },
         i = {
            ali[(nseq1+1):nrow(ali), n] <- pdbs2$ali[, k]
            xyz[(nseq1+1):nrow(ali), atom2xyz(n)] <- pdbs2$xyz[, atom2xyz(k)]
            k = k-1 
         },
         d = {
            ali[1:nseq1, n] <- pdbs1$ali[, j]
            xyz[1:nseq1, atom2xyz(n)] <- pdbs1$xyz[, atom2xyz(j)]
            j = j-1 
         })
      n = n-1
   }
   ali <- ali[, !apply(ali, 2, function(x) all(is.gap(x)))]  # trim common gaps
   xyz <- xyz[, !apply(xyz, 2, function(x) all(is.na(x)))]  # trim common gaps
   rownames(ali) <- c(rownames(pdbs1$ali), rownames(pdbs2$ali))
   rownames(xyz) <- c(rownames(pdbs1$xyz), rownames(pdbs2$xyz))
   aln <- list(id=c(pdbs1$id, pdbs2$id), ali=ali, xyz=xyz)
   if(!ret.score)
      return( aln )
   else
      return( list(aln=aln, score=mat[nrow(score.mat)+1, ncol(score.mat)+1]) )
}
