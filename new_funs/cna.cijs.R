# ensemble CNA calculation optimized with multicore
cna.cijs <- function(cijs, ..., ncore = NULL) {
   ncore <- setup.ncore(ncore)
   
   if("all.dccm" %in% names(cijs)) cijs <- cijs$all.dccm 
   if(is.matrix(cijs)) {
      net <- cna(cijs, ...)
   } else {
      if(is.array(cijs) && length(dim(cijs)==3))
         cijs <- do.call("c", apply(cijs, 3, list))
      if(is.list(cijs)) {
         net <- mclapply(cijs, cna, ...)
      } else {
         warning("cijs should be matrix, array(dim=3), or list")
         net <- NULL 
      }
   }
   return(net)
}
