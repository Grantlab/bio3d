layout.pdb <- function(pdb, membership, renumber=TRUE, k=3){

  ## Return the coordinate centers of input communitys
  ##  as defined in 'membership vector' using Calpha's in 'pdb'.
  ##  If k=3 the xyz geometric centers are returned. if k<3 then
  ##  multidimensional scaling is used for k space ordination.
  ##
  ##   co2 <- layout.pdb(pdb, net$raw.communities$membership, k=2)
  ##   co2 <- layout.pdb(pdb, net, k=2)
  ##   plot.cna(net, layout=co2)

  if( class(pdb) != "pdb" ) {
    stop("Input 'pdb' is not of class 'pdb' as obtained from 'read.pdb()'")
  }
  if(!k %in% c(1,2,3)) {
    stop("Input 'k' should have a value of 3, 2 or 1")
  }
  if((class(membership) == "cna") || is.list(membership)) {
    membership <- membership$communities$membership
  }

  ## Renumber 'pdb' to match membership resno indices
  if(renumber) {
    pdb <- convert.pdb(pdb, type="pdb", renumber=TRUE, rm.h=FALSE) ## verbose=FALSE
  }

  ##-- Check if the number of number of residues in 'pdb' equals
  ##   the length of 'membership' vector
  notprotein.inds <- atom.select(pdb, "notprotein", verbose=FALSE)

  if(length(notprotein.inds$atom)>0){
    num.res <- length(pdb$atom[pdb$calpha,"resno"]) + length(unique(pdb$atom[notprotein.inds$atom,6]))
  }
  if(length(notprotein.inds$atom)==0){
    num.res <- length(pdb$atom[pdb$calpha,"resno"])
  }
  
  ##-- Calculate the geometric center of each community
  n <- unique(membership[!is.na(membership)])
  
  cent <- matrix(NA, nrow=length(n), ncol=3)
  a <- 1
  for(i in n){
    inds <- atom.select(pdb, resno=which(membership==i),
                        elety="CA", verbose=FALSE)

    cent[a,] <- apply( matrix(pdb$xyz[inds$xyz], nrow=3), 1, mean)
    a <- a + 1 
  }

  if(k != 3) {
    ##-- Multidimensional dcaling for 2D or 1D projection
    ##   note. dist(centers) and dist.xyz(centers) give same answer
    cent <- cmdscale(dist(cent),k=k)
  }
  return(cent)
}
