"get.pdb" <-
  function (ids, path = "./", URLonly = FALSE, overwrite = FALSE ) 
{
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
    for (k in 1:length(pdb.files)) {
      
      if (!file.exists(put.files[k]) | overwrite ) {
        rtn[k] <- download.file(get.files[k], put.files[k])
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
