local.fitting <- function(trj, pdb, radius=9, ncore=1){

  if(!is.numeric(trj)){
    stop("Trj object must be a numeric matrix such as obtained from read.ncdf")
  }

  if(sum(class(pdb) %in% "pdb")==0){
    stop("Pdb must be a pdb object such as obtained from read.pdb")
  }

  if(!is.numeric(radius)){
    stop("Radius must a numeric value")
  }

  ca.inds <- atom.select(pdb, elety="CA")
  rmsd.matrix <- matrix(NA, nrow=(dim(trj)[1]-1), ncol=length(ca.inds$atom))
  
  for(i in 2:dim(trj)[1]){
    print(paste0("frame analysed ",i))
    ## calculation of CA distances
    CA.dist <- cmap(trj[i,ca.inds$xyz], dcut=radius, scut=0, ncore=ncore,mask.lower=FALSE)
  
    for(j in 1:length(ca.inds$atom)){
      #print(j)
      neighbors <- which(CA.dist[,j] == 1)
      neighbors.inds <- atom2xyz(neighbors)
      
      coords.fitted <- fit.xyz(trj[i-1,], trj[i,], fixed.inds=ca.inds$xyz[neighbors.inds], mobile.inds=ca.inds$xyz[neighbors.inds])

      rmsd.matrix[i-1,j] <- rmsd(coords.fitted, trj[i-1,], a.inds=ca.inds$xyz[neighbors.inds], b.inds=ca.inds$xyz[neighbors.inds])
    }
  }

  correlation.matrix <- matrix(NA, nrow=length(ca.inds$atom), ncol=length(ca.inds$atom))

  for(i in 1:(length(ca.inds$atom)-1)){
    for(j in (i+1):length(ca.inds$atom)){
      correlation.matrix[i,j] <- cor(rmsd.matrix[,i], rmsd.matrix[,j], method="pearson")
    }
  }

  output <- list("correlation.matrix"=correlation.matrix, "rmsd.matrix"=rmsd.matrix)
  
  return(output)
}
      
  
