dccm.local.xyz <- function(trj, pdb, radius=9, basin.grp=NULL, write.average=FALSE, ncore=1){

  ## Check for multiple cores
  ncore = setup.ncore(ncore)

  if(ncore > 1) {
    require(parallel)
    mcparallel <- get("mcparallel", envir = getNamespace("parallel"))
  }

  if((!is.numeric(trj)) || (!is.matrix(trj))){
    stop("trj object must be a numeric atomic coordinates matrix such as obtained from read.ncdf")
  }

  if(sum(class(pdb) %in% "pdb")==0){
    stop("pdb must be a pdb object such as obtained from read.pdb")
  }

  if(!is.numeric(radius)){
    stop("radius must a numeric value")
  }
  
  
  ## Internal function to do local fit and rmsd calculations
  fragment.rmsd <- function(x, ref, pdb, ca.inds, radius, CA.map=NULL) {
  
      ## calculation of CA distances
      if(is.null(CA.map)){
        CA.map <- cmap(ref[ca.inds$xyz], dcut=radius, scut=0, mask.lower=FALSE)
      }
      
      sapply(1:length(ca.inds$atom), function(j) {
     
        neighbors <- which(CA.map[,j] == 1)
        resno = pdb$atom[ca.inds$atom, "resno"][neighbors]
        neighbors.inds <- atom.select(pdb, resno=resno, string="noh", verbose=FALSE)
      
        coords.fitted <- fit.xyz(ref, x, fixed.inds=neighbors.inds$xyz, mobile.inds=neighbors.inds$xyz)

        rmsd(coords.fitted, ref, a.inds=neighbors.inds$xyz, b.inds=neighbors.inds$xyz)

      } )
  }
   
  ca.inds <- atom.select(pdb, elety="CA")
 
   ## Fitting and RMSD on previuos frame
   if(is.null(basin.grp)){

     ## Initialize progress bar
     pb <- txtProgressBar(min=1, max=dim(trj)[1]-1, style=3)

     if(ncore > 1) {   # Parallel

       # For progress bar
       fpb <- fifo(tempfile(), open = "w+b", blocking = T)

       # spawn a child process for message printing
       child <- mcparallel({
         progress <- 0.0
         while(progress < (dim(trj)[1]-1) && !isIncomplete(fpb)) {
           msg <- readBin(fpb, "double")
           progress <- progress + as.numeric(msg)
           setTxtProgressBar(pb, progress)
         }
       } )
      ###################

       rmsd.list <- mclapply(2:dim(trj)[1], function(i) {
         writeBin(1, fpb)
         fragment.rmsd(trj[i, ], trj[i-1, ], pdb, ca.inds, radius)
         }) 
       rmsd.matrix = do.call(rbind, rmsd.list)

       close(fpb)
     }  
     else {

       rmsd.matrix <- t( sapply(2:dim(trj)[1], function(i) {
         setTxtProgressBar(pb, i-1)
         fragment.rmsd(trj[i, ], trj[i-1, ], pdb, ca.inds, radius)
       }) )
     }
   }
   
   if(!is.null(basin.grp)){

     if(!is.numeric(basin.grp)){
       stop("basin.grp must a numeric value")
     }

     if(length(basin.grp) != dim(trj)[1]){
       stop("basin.grp must be of the same length as the number of rows in trj")
     }
     
     ## average structures calculation
     average <- matrix(NA, nrow=length(unique(basin.grp)), ncol=dim(trj)[2])
     for(h in 1:length(unique(basin.grp))){
       average[h,] <- apply(trj[c(which(basin.grp==h)),], 2, mean)

       if(sum(write.average)>0){
         write.pdb(pdb, xyz=average[h,], file=paste0("average.basin.", h, ".pdb"))
       }

     }
     
     ## Initialize progress bar
     pb <- txtProgressBar(min=1, max=dim(trj)[1]-1, style=3)

     if(ncore > 1) {   # Parallel

       # For progress bar
       fpb <- fifo(tempfile(), open = "w+b", blocking = T)

       # spawn a child process for message printing
       child <- mcparallel({
         progress <- 0.0
         while(progress < (dim(trj)[1]-1) && !isIncomplete(fpb)) {
           msg <- readBin(fpb, "double")
           progress <- progress + as.numeric(msg)
           setTxtProgressBar(pb, progress)
         }
       } )
      ###################

       rmsd.list.tmp <- list()
       for(h in 1:dim(average)[1]){
         
         CA.map <- cmap(average[h,ca.inds$xyz], dcut=radius, scut=0, mask.lower=FALSE, ncore=ncore)

         rmsd.list.tmp <- c(rmsd.list.tmp, mclapply(min(which(basin.grp==h)):max(which(basin.grp==h)), function(i) {
           writeBin(1, fpb)
           fragment.rmsd(trj[i, ], average[h,], pdb, ca.inds, radius)
           })) 
       }
       
       rmsd.matrix = do.call(rbind, rmsd.list.tmp)

       close(fpb)
     }  
     else {

       rmsd.matrix <- NULL
       for(h in 1:dim(average)[1]){
         CA.map <- cmap(average[h,ca.inds$xyz], dcut=radius, scut=0, mask.lower=FALSE)
         
         rmsd.matrix <- rbind(rmsd.matrix, t( sapply(min(which(basin.grp==h)):max(which(basin.grp==h)), function(i) {
           setTxtProgressBar(pb, i-1)
           fragment.rmsd(trj[i, ], average[h,], pdb, ca.inds, radius, CA.map=CA.map)
         }) ))
       }

     }
   }
   
   close(pb)
   
   correlation.matrix <- cor(rmsd.matrix)

   output <- list("correlation.matrix"=correlation.matrix, "rmsd.matrix"=rmsd.matrix)
  
  return(output)
}
      
  

# basin.grp <- c(rep(1,495), rep(2,565), rep(3,416), rep(4,392), rep(5,132))
