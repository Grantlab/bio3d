"get.pdb" <-
  function (ids, path = ".", URLonly = FALSE, overwrite = FALSE, verbose = TRUE, ncore = 1 ) 
{
    # Parallelized by multicore package (Tue Oct 15 15:23:36 EDT 2013)
    if(ncore > 1) {
       oops <- require(multicore)
       if(!oops)
          stop("Please install the multicore package from CRAN")
       if(ncore > 8) {
          # To avoid too frequent access to PDB server
          warning("Maximum ncore (=8) exceeds")
          ncore = 8
       }
       options(cores = ncore)
    }

    if (any(nchar(ids) != 4)) {
        warning("ids should be standard 4 character PDB formart: trying first 4 char...")
        ids <- substr(basename(ids), 1, 4)
    }
    ids <- unique(ids)
    pdb.files <- paste(ids, ".pdb", sep = "")
    get.files <- file.path("http://www.rcsb.org/pdb/files", pdb.files)
    if (URLonly) 
        return(get.files)
    put.files <- file.path(path, pdb.files)

    if(!file.exists(path))
       dir.create(path)
    rtn <- rep(NA, length(pdb.files))
    if(ncore>1) {
       rtn <- unlist(mclapply(1:length(pdb.files), function(k) {
          if (!file.exists(put.files[k]) | overwrite ) {
            rtn <- download.file(get.files[k], put.files[k], quiet = !verbose)
          }
          else {
            rtn <- put.files[k]
            warning(paste(put.files[k], " exists. Skipping download"))
          }
          return(rtn)
       }))
       readChildren()
    } else {
       for (k in 1:length(pdb.files)) {
         
         if (!file.exists(put.files[k]) | overwrite ) {
           rtn[k] <- download.file(get.files[k], put.files[k], quiet = !verbose)
         }
         else {
           rtn[k] <- put.files[k]
           warning(paste(put.files[k], " exists. Skipping download"))
         }
         
       }
    }
    names(rtn) <- file.path(path, paste(ids, ".pdb", sep = ""))
    if (any(rtn == 1)) {
        warning("Some files could not be downloaded, check returned value")
        return(rtn)
    }
    else {
        return(names(rtn))
    }
}
