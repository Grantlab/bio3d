multiple.sims.network <- function(cij=cij,
                                  cmap=cmap,
                                  vmd.out=NULL){

  net.no.cmap.btwn.1 <- cna(cij)
  net.no.cmap.btwn.2 <- cna(cij, cutoff.cij=0.6)
  net.cmap.btwn.1 <- cna(cij, cutoff.cij=0.2, cm=cmap)
  net.cmap.btwn.2 <- cna(cij, cm=cmap)
  net.cmap.btwn.3 <- cna(cij, cutoff.cij=0.6, cm=cmap)
  net.no.cmap.walk.1 <- cna(cij, cluster.method="walk")
  net.no.cmap.walk.2 <- cna(cij, cutoff.cij=0.6, cluster.method="walk")
  net.cmap.walk.1 <- cna(cij, cutoff.cij=0, cluster.method="walk", cm=cmap)
  net.cmap.walk.2 <- cna(cij, cluster.method="walk", cm=cmap)
  net.cmap.walk.3 <- cna(cij, cutoff.cij=0.6, cluster.method="walk", cm=cmap)
  net.no.cmap.greed.1 <- cna(cij, cluster.method="greed")
  net.no.cmap.greed.2 <- cna(cij, cutoff.cij=0.6, cluster.method="greed")
  net.cmap.greed.1 <- cna(cij, cluster.method="greed", cm=cmap)
  net.cmap.greed.2 <- cna(cij, cutoff.cij=0.6, cluster.method="greed", cm=cmap)

  cij.vmd <- cij * cmap
  cij.vmd[abs(cij.vmd)<0.4] <- 0
  diag(cij.vmd) <- 0
  cij.vmd[cij.vmd>0] <- -log(cij.vmd[cij.vmd>0])
  net.cmap.vmd <- vmd.comms(vmd.out,adjmatrix=cij.vmd)

  output <- list("net.no.cmap.btwn.1"=net.no.cmap.btwn.1,
                 "net.no.cmap.btwn.2"=net.no.cmap.btwn.2,
                 "net.cmap.btwn.1"=net.cmap.btwn.1,
                 "net.cmap.btwn.2"=net.cmap.btwn.2,
                 "net.cmap.btwn.3"=net.cmap.btwn.3,
                 "net.no.cmap.walk.1"=net.no.cmap.walk.1,
                 "net.no.cmap.walk.2"=net.no.cmap.walk.2,
                 "net.cmap.walk.1"=net.cmap.walk.1,
                 "net.cmap.walk.2"=net.cmap.walk.2,
                 "net.cmap.walk.3"=net.cmap.walk.3,
                 "net.no.cmap.greed.1"=net.no.cmap.greed.1,
                 "net.no.cmap.greed.2"=net.no.cmap.greed.2,
                 "net.cmap.greed.1"=net.cmap.greed.1,
                 "net.cmap.greed.2"=net.cmap.greed.2,
                 "net.cmap.vmd"=net.cmap.vmd)
                 

  return(output)
}       
