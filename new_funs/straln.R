# Multiple Structure Alignment Based on
# Dynamic Programming and Tree-guided Progressive Alignment
# Currently fitting is not supported... (Assume pre-fitted structures with some reference)

merge.aln <- function(aln1, aln2) {
   aln <- NULL
   aln$id <- c(aln1$id, aln2$id)
   aln$ali <- rbind(aln1$ali, aln2$ali)
   if(!is.null(aln1$xyz))
      aln$xyz <- rbind(aln1$xyz, aln2$xyz)
   return( aln )
}

sub.aln <- function(aln, inds=1:nrow(aln$ali), pos=1:ncol(aln$ali)) {
   if(length(inds)==0 || any(inds<=0) ||
      any(is.na(inds)) ) stop("Error: invalid inds value")
   naln <- NULL
   naln$id <- aln$id[inds]
   naln$ali <- rbind(NULL, aln$ali[inds, pos])
   naln$ali <- rbind(NULL, naln$ali[, !apply(naln$ali, 2, function(x) all(is.gap(x)))])  # trim common gaps
   if(!is.null(aln$xyz)) {
      naln$xyz <- rbind(NULL, aln$xyz[inds, atom2xyz(pos)])
      naln$xyz <- rbind(NULL, naln$xyz[, !apply(naln$xyz, 2, function(x) all(is.na(x)))])  # trim common gaps
   }
   names <- rownames(aln$ali)
   if(!is.null(names)) rownames(naln$ali) <- names[inds]
   names <- rownames(aln$xyz)
   if(!is.null(names)) rownames(naln$xyz) <- names[inds]
   return( naln )
}

straln <- function(pdbs, prefix = "", outfile = "straln.fa", ret.hc = FALSE, ncore=NULL, ...) {
   ncore <- setup.ncore(ncore)

   dots <- list(...)
   if("gap.penalty" %in% names(dots)) gap.penalty = dots$gap.penalty
   else gap.penalty = -9

#   gap.penalty <- - gap.penalty^2  #MSF

   mylapply <- lapply
   if(ncore > 1) mylapply <- mclapply
  
   # Full pairwise compairsion 
   np <- pairwise(nrow(pdbs$ali))
   dist.list <- mylapply(1:nrow(np), function(i) {
      ii = np[i, 1]
      jj = np[i, 2]
      pdbs1 <- sub.aln(pdbs, ii) 
      pdbs2 <- sub.aln(pdbs, jj) 
      taln <- straln.pair(pdbs1, pdbs2, ret.score=TRUE, ncore=1, ...)
      gaps <- gap.inspect(taln$aln$ali)
      rmsd <- sqrt( - (taln$score - gap.penalty*length(gaps$t.inds))/length(gaps$f.inds) )
#      rmsd <- sqrt( - taln$score/ncol(taln$aln$ali) )
      return( rmsd )
   })
   dd <- unlist(dist.list)
   distmat <- matrix(0, nrow(pdbs$ali), nrow(pdbs$ali))
   distmat[lower.tri(distmat)] <- dd
   distmat[upper.tri(distmat)] <- t(distmat)[upper.tri(distmat)]

   #phy <- nj(as.dist(distmat))
   hc <- hclust(as.dist(distmat))  # Guide tree

   # Progressive alignment
   aln <- NULL
   for(i in 1:nrow(hc$merge)) {
      ii = hc$merge[i, 1]
      jj = hc$merge[i, 2]
      if(ii < 0) pdbs1 <- sub.aln(pdbs, -ii)
      else pdbs1 <- aln[[ii]]
      if(jj < 0) pdbs2 <- sub.aln(pdbs, -jj)
      else pdbs2 <- aln[[jj]]

      aln <- c(aln, list(straln.pair(pdbs1, pdbs2, ncore=ncore, ...)))
   }
   aln <- aln[[nrow(hc$merge)]]
   inds <- match(pdbs$id, aln$id)
   aln <- list(id=aln$id[inds], ali=aln$ali[inds, ])
   if(!is.null(outfile)) write.fasta(aln, file = outfile)
   if(ret.hc) return( list(aln=aln, hc=hc, distmat=distmat) )
   else return( aln )
}
