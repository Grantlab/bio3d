`cmap` <-
function(xyz, grpby=NULL, dcut=4, scut=3, mask.lower = TRUE) {

  ## Distance matrix (all-atom)
  dmat <- dm.xyz( xyz, grpby, scut, mask.lower = mask.lower)
  ## Contact map
  return(matrix(as.numeric(dmat < dcut),
                ncol = ncol(dmat),
                nrow = nrow(dmat)))
}

