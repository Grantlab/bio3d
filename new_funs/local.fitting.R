local.fitting <- function(trj, pdb, radius=9, basin.grp=NULL, ncore=1){

   ## Check for multiple cores
   ncore = setup.ncore(ncore)

   if(ncore > 1) {
     require(parallel)
     mcparallel <- get("mcparallel", envir = getNamespace("parallel"))
  }

  if(!is.numeric(trj)){
    stop("Trj object must be a numeric matrix such as obtained from read.ncdf")
  }

  if(sum(class(pdb) %in% "pdb")==0){
    stop("Pdb must be a pdb object such as obtained from read.pdb")
  }

  if(!is.numeric(radius)){
    stop("Radius must a numeric value")
  }
   
  ## Internal function to do local fit and rmsd calculations
  fragment.rmsd <- function(x, ref, pdb, ca.inds, radius) {
  
      ## calculation of CA distances
      CA.map <- cmap(x[ca.inds$xyz], dcut=radius, scut=0, mask.lower=FALSE)
      
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
  close(pb)
   
  correlation.matrix <- cor(rmsd.matrix)
    

  output <- list("correlation.matrix"=correlation.matrix, "rmsd.matrix"=rmsd.matrix)
  
  return(output)
}
      
  
