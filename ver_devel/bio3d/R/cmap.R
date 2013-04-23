`cmap` <-
function(xyz, grpby=NULL, dcut=4, scut=3, pcut=1, mask.lower = TRUE,
         ncore=1, nseg.scale=1) {

  # Parallelized by multicore package (Mon Apr 22 16:32:19 EDT 2013)
  if(ncore > 1) {
     oops <- require(multicore)
     if(!oops)
        stop("Please install the multicore package from CRAN")

     options(cores = ncore)

     # Issue of serialization problem
     # Maximal number of cells of a double-precision matrix
     # that each core can serialize: (2^31-1-61)/8
     R_NCELL_LIMIT_CORE = 2.68435448e8
     R_NCELL_LIMIT = ncore * R_NCELL_LIMIT_CORE

     if(nseg.scale < 1) {
        warning("nseg.scale should be 1 or a larger integer\n")
        nseg.scale=1
     }
  }

  if(is.matrix(xyz)) {
     if(is.null(grpby)) {
        nres <- ncol(xyz)/3
     } else {
        inds <- bounds(grpby, dup.inds = TRUE)
        nres <- nrow(inds)
     }
     if(ncore > 1) {
        ni = nrow(xyz)
        RLIMIT = R_NCELL_LIMIT
        nDataSeg = floor((ni-1)/RLIMIT)+1
        nDataSeg = floor(nDataSeg * nseg.scale)
        lenSeg = floor(ni/nDataSeg)
        cmap.list <- NULL
        for(i in 1:nDataSeg) {
           istart = (i-1)*lenSeg + 1
           iend = if(i<nDataSeg) i*lenSeg else ni
           cmap.list <- c(cmap.list, mclapply(istart:iend, function(j) {
               dmat <- dm.xyz(xyz[j,], grpby, scut, mask.lower=TRUE)
               return(as.numeric(dmat[!lower.tri(dmat)] < dcut))
           }) )
        }
        readChildren()
     } else {
        cmap.list <- lapply(1:nrow(xyz), function(j) {
            dmat <- dm.xyz(xyz[j,], grpby, scut, mask.lower=TRUE)
            return(as.numeric(dmat[!lower.tri(dmat)] < dcut))
        }) 
     }
     cmap.t <- rowMeans(do.call(cbind, cmap.list))
     cmap.t <- as.numeric(cmap.t >= pcut )
     cont.map <- matrix(NA, nrow=nres, ncol=nres)
     cont.map[!lower.tri(cont.map)] <- cmap.t
     if(!mask.lower) 
         cont.map[lower.tri(cont.map)] <- t(cont.map)[lower.tri(cont.map)]
  } else {
     ## Distance matrix (all-atom)
     dmat <- dm.xyz( xyz, grpby, scut, mask.lower = mask.lower)
     ## Contact map
     return(matrix(as.numeric(dmat < dcut),
                ncol = ncol(dmat),
                nrow = nrow(dmat)))
  }
  return (cont.map)
}

