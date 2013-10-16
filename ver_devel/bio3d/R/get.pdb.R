"get.pdb" <-
  function (ids, path = ".", URLonly = FALSE, overwrite = FALSE, gzip = FALSE ) 
{
    if (any(nchar(ids) != 4)) {
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
    
    for (k in 1:length(pdb.files)) {
      if (!file.exists(sub(".gz$", "", put.files[k])) | overwrite ) {
        rtn[k] <- download.file(get.files[k], put.files[k])
        if(gzip) {
          cmd <- paste("gunzip -f", put.files[k])
          os1 <- .Platform$OS.type
          if (os1 == "windows")
            shell(shQuote(cmd))
          else 
            system(cmd)
        }
      }
      else {
        rtn[k] <- put.files[k]
        warning(paste(put.files[k], " exists. Skipping download"))
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
