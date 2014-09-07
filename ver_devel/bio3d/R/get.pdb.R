"get.pdb" <-
  function (ids, path = ".", URLonly = FALSE, overwrite = FALSE, gzip = FALSE, verbose = TRUE, ncore = 1 ) 
{
    if(.Platform$OS.type=="windows")
      gzip <- FALSE
    
    # Parallelized by parallel package (Tue Oct 15 15:23:36 EDT 2013)
    ncore <- setup.ncore(ncore)
    if(ncore > 4) {
       # To avoid too frequent access to PDB server
       warning("Exceed maximum ncore (=4) to access PDB server. Use ncore=4")
       ncore <- setup.ncore(ncore = 4)
    }

    if(inherits(ids, "blast")) ids = ids$pdb.id

    if (any(nchar(ids) < 4)) stop("ids should be standard 4 character PDB-IDs")
    if (any(nchar(ids) > 4)) {
        warning("ids should be standard 4 character PDB-IDs: trying first 4 characters...")
        ids <- substr(basename(ids), 1, 4)
    }
    ids <- unique(ids)
    pdb.files <- paste(ids, ".pdb", ifelse(gzip, ".gz", ""), sep = "")
    get.files <- file.path("http://www.rcsb.org/pdb/files", pdb.files)
    if (URLonly) 
        return(get.files)
    put.files <- file.path(path, pdb.files)

    if(!file.exists(path))
       dir.create(path)
    rtn <- rep(NA, length(pdb.files))

    if(ncore > 1) {
       rtn <- unlist(mclapply(1:length(pdb.files), function(k) {
          if (!file.exists(sub(".gz$", "", put.files[k])) | overwrite ) {
            rtn <- try(download.file(get.files[k], put.files[k], quiet = !verbose), silent = TRUE)
            if(inherits(rtn, "try-error")) {
               rtn <- 1
               file.remove(put.files[k])
            } else if(gzip) {
               cmd <- paste("gunzip -f", put.files[k])
               system(cmd)
            }
          }
          else {
            rtn <- put.files[k]
            warning(paste(put.files[k], " exists. Skipping download"))
          }
          return(rtn)
       }))
    } else {
       for (k in 1:length(pdb.files)) {
         if (!file.exists(sub(".gz$", "", put.files[k])) | overwrite ) {
           rt <- try(download.file(get.files[k], put.files[k], quiet = !verbose), silent=TRUE)
           rtn[k] <- rt
           if(inherits(rt, "try-error")) {
              rtn[k] <- 1
              file.remove(put.files[k])
           } else if(gzip) {
              cmd <- paste("gunzip -f", put.files[k])
              system(cmd)
           }
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
