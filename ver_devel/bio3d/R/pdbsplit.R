`pdbsplit` <-
function(pdb.files, ids=NULL, path="split_chain", multi=FALSE) {
  out <- c()
  toread <- file.exists(pdb.files)
  toread[substr(pdb.files, 1, 4) == "http"] <- TRUE
  if (all(!toread)) 
    stop("No corresponding PDB files found")
  if (any(!toread)) {
    warning(paste("Missing files:\n\t",
                  paste(pdb.files[!toread], 
                        collapse = "\n\t"), sep = ""))
    
    pdb.files <- pdb.files[toread]
  }
  
  if(!file.exists(path)) 
     dir.create(path)
  for (i in 1:length(pdb.files)) {
    pdb <- read.pdb(pdb.files[i], multi=multi, het2atom = TRUE, maxlines = -1)
    chains <- unique(pdb$atom[, "chain"])
    if (length(chains) > 0) {
      for (j in 1:length(chains)) {
        if (!is.na(chains[j])) {
          new.pdb <- NULL
          
          sel <- atom.select(pdb, paste("//", chains[j], "/////"))
          new.pdb <- trim.pdb(pdb, sel)

          ## Multi-model records
          if (!is.null(pdb$xyz.models)) {

            for ( k in 1:nrow(pdb$xyz.models) ) {
              
              str.len <- nchar(nrow(pdb$xyz.models))
              new.name <- paste(substr(basename(pdb.files[i]), 
                                       1, 4), "_", chains[j], ".",
                                formatC(k, width=str.len, format="d", flag="0"),
                                ".pdb", sep = "")
              new.name <- file.path(path, new.name)
              
              xyz <- pdb$xyz.models[k, sel$xyz]
              write.pdb(new.pdb, file = new.name, xyz=xyz)
              out <- c(out, new.name)
            }
          }
          else {
            new.name <- paste(substr(basename(pdb.files[i]), 
                                     1, 4), "_", chains[j], ".pdb", sep = "")
            new.name <- file.path(path, new.name)
            write.pdb(new.pdb, file = new.name)
            out <- c(out, new.name)
          }
        }
      }
    }
  }

  if(!is.null(ids))
    out <- out[ unlist(lapply(ids, grep, out)) ]
  return(out)
}
