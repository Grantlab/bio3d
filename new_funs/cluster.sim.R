## cluster.sim()
##
## Date:     Sun Sep 15 19:54:49 EDT 2013
##
## Purpose:  Computation of Similarity Measures for Clustering Results 
## Notes:    Adapted from cindex() function of Thomas Girke
## Utility:  Similarity coefficients for two clustering results C and K. 
## The current set of methods includes the Jaccard Index, the Rand 
## Index and the Adjusted Rand Index.
## 
## Definition of Indices:
##   - Jaccard Index: J(C,K) = a/(a+b+c) 
##   - Rand Index: R(C,K) = (a+d)/(a+b+c+d) 
##   - Adjusted Rand Index: AR(C,K) = (2*(a*d-c*b)) / ((a+b)*(b+d)+(a+c)*(c+d)) 
##      Where:
##      a = number of pairs that cluster together in both C and K
##      b = number of pairs that cluster together in C, but not K
##      c = number of pairs that cluster together in K, but not C
##      d = number of pairs that are not joined into clusters in K and C
##
## Useage: 
##	cindex(grp1, grp22, self=FALSE, minSZ=1, method="jaccard")
##               # grp1/grp2: named vectors where the cluster IDs are values and the item labels are names 
##               # method: supported similarity methods are "jaccard", "rand" and "arand" 
##               # self: FALSE to ignore or TRUE to include clusters with single items
##               # minSZ: selection of a minimum cluster size
##
## More detailed instructions are available on this page:
##     http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/R_BioCondManual.html#clustering_jaccard


cluster.sim <- function(grp1, grp2, self=FALSE, minSZ=1, method="jaccard") {
  ## Compute Jaccard and Rand Indices for two clustering results
  ## (grp1 and grp2 are named vectors of cluster IDs (values) of residues (names).

  if(is.null(names(grp1)) || is.null(grp2) ) {
    stop("Input cluster membership verctors sould have variables/residues as names attribute")
  }
  
  ## (A) Conversions to assure proper data structure 
  ## Allow selection of a minimum cluster size. 
  if(minSZ > 1) {
    CLSZ <- table(grp1)
    grp1 <- grp1[grp1 %in% names(CLSZ[CLSZ >= minSZ])]
    CLSZ <- table(grp2)
    grp2 <- grp2[grp2 %in% names(CLSZ[CLSZ >= minSZ])]
  }
  
  ## Remove non-common labels in grp1 and grp2 
  grp1 <- grp1[names(grp1) %in% names(grp2)]
  grp2 <- grp2[names(grp2) %in% names(grp1)]
  
  ## (B) Compute matrix that provides the common pairs among the two clustering results
  cpairMA <- sapply(unique(grp1), function(x)  sapply(unique(grp2), function(z) sum(names(grp1[grp1==x]) %in% names(grp2[grp2==z]))))
  
  colnames(cpairMA) <- unique(grp1)
  rownames(cpairMA) <- unique(grp2)
  olMA <- cpairMA # Copy for export.
  ### To debug result, use: sum(names(grp1[grp1==2]) %in% names(grp2[grp2==28]))

  
  if(self==TRUE) {
    ## Plus cpairMA/2 includes self comparisons
    cpairMA <- cpairMA^2/2 + cpairMA/2 
  }
  if(self==FALSE) {
    ## Minus cpairMA/2 excludes self comparisons
    cpairMA <- cpairMA^2/2 - cpairMA/2 
  }

  
  ## (C) Compute Similarity Indices
  if(self==TRUE) { 
    NpairsCL1 <- tapply(names(grp1), grp1, length)
    NpairsCL1 <- NpairsCL1^2/2 + NpairsCL1/2
    
    NpairsCL2 <- tapply(names(grp2), grp2, length)
    NpairsCL2 <- NpairsCL2^2/2 + NpairsCL2/2
  }
  if(self==FALSE) { 
    NpairsCL1 <- tapply(names(grp1), grp1, length)
    NpairsCL1 <- NpairsCL1^2/2 - NpairsCL1/2

    NpairsCL2 <- tapply(names(grp2), grp2, length)
    NpairsCL2 <- NpairsCL2^2/2 - NpairsCL2/2
  }
  
  ja <- sum(cpairMA)
  jb <- sum(NpairsCL1[colnames(cpairMA)] - colSums(cpairMA))
  jc <- sum(NpairsCL2[rownames(cpairMA)] - rowSums(cpairMA))

  ## Return results as list
  if(method=="jaccard") {
    jindex <- ja/(ja+jb+jc)
    
    return(list(intersects=olMA,
                variables=unlist(list(a=ja, b=jb, c=jc)),
                Jaccard_Index=jindex))
  }
  if(method=="rand") {
    ## grp1 and grp2 contain same items (see above)
    Nitems <- length(names(grp1))
    
    if(self==TRUE) {
      jd <- (Nitems^2/2 + Nitems/2) - ja - jb -jc
    }
    if(self==FALSE) {
      jd <- (Nitems^2/2 - Nitems/2) - ja - jb -jc
    }
    rindex <- (ja+jd)/(ja+jb+jc+jd)	

    return(list(intersects=olMA,
                variables=unlist(list(a=ja, b=jb, c=jc, d=jd)),
                Rand_Index=rindex))
  }
  if(method=="arand") {
    Nitems <- length(names(grp1))
    if(self==TRUE) {
      jd <- (Nitems^2/2 + Nitems/2) - ja - jb -jc
    }
    if(self==FALSE) {
      jd <- (Nitems^2/2 - Nitems/2) - ja - jb -jc
    }
    arindex <- (2*(ja*jd - jc*jb)) / ((ja+jb)*(jb+jd)+(ja+jc)*(jc+jd))	

    return(list(intersects=olMA,
                variables=unlist(list(a=ja, b=jb, c=jc, d=jd)),
                Adjusted_Rand_Index=arindex))
  }
}
