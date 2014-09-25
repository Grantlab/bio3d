prune.couplings <- function(x, sse, neigh=2){

  if(class(sse) != "sse"){
    stop("sse object must be have class 'sse', such as obtained from 'dssp'")
  }

  if(!is.numeric(x)){
    stop("x object must be numeric")
  }

  if(dim(x)[1] != dim(x)[2]){
    stop("x object must be a square matrix")
  }
  
  ## pruning alpha helix couplings
  num.of.helices <- length(sse$helix$start)

  for(i in 1:num.of.helices){
    helix.residues <- c(sse$helix$start[i]:sse$helix$end[i])

    num.of.helix.residues <- length(helix.residues)

    if(num.of.helix.residues > neigh){
      for(j in c(helix.residues[1]:helix.residues[num.of.helix.residues-neigh])){
        for(h in c((j+neigh):helix.residues[num.of.helix.residues])){
          x[j,h] <- 0
        }
      }
    }
  }

  ## pruning beta sheet couplings
  num.of.sheet <- length(sse$sheet$start)

  for(i in 1:num.of.sheet){
    sheet.residues <- c(sse$sheet$start[i]:sse$sheet$end[i])

    num.of.sheet.residues <- length(sheet.residues)

    if(num.of.sheet.residues > neigh){
      for(j in c(sheet.residues[1]:sheet.residues[num.of.sheet.residues-neigh])){
        for(h in c((j+neigh):sheet.residues[num.of.sheet.residues])){
          x[j,h] <- 0
        }
      }
    }
  }

  temp <-  t(x)
  x[lower.tri(x)] <- temp[lower.tri(temp)]

  return(x)
}
