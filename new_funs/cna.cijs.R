# ensemble CNA calculation optimized with multicore
cna.cijs <- function(cijs, ..., ncore = NULL) {
   require(parallel)
   if(is.null(ncore)) ncore = detectCores()
   options(mc.cores=ncore)
   
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
