`pdbsplit` <-
function(pdb.files, ids=NULL, path="split_chain", ...) {
  out <- c(); unused <- c()
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
    pdb <- read.pdb(pdb.files[i], ...)
    chains <- unique(pdb$atom[, "chain"])
    
    if(!is.null(ids)) {
      ids <- unique(ids)
      
      ## Match 'ids' with 'pdbId_chainId' combinations
      tmp.names <- paste(substr(basename(pdb.files[i]), 
                                1, 4), "_", chains, sep = "")

      tmp.inds <- unique(unlist(lapply(ids, grep, tmp.names)))
      if(length(tmp.inds)==0) {
        ## Skip pdb file if no match were found
        unused <- c(unused, substr(basename(pdb.files[i]), 1, 4))
        chains <- c()
      }
      else {
        chains <- chains[tmp.inds]
      }
    }
    
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
  
  if(!is.null(ids)) {
    ids.used <- NULL; nonmatch <- NULL
    if(length(out)>0) {
      ids.used <- sub(".pdb$", "", basename(out))
      tmp.fun <- function(x, y) { ifelse(length(grep(x,y))>0, TRUE, FALSE) }
      tmp.inds <- unlist(lapply(ids, tmp.fun, ids.used))
      nonmatch <- ids[!tmp.inds]
    }

    ## Elements of 'pdb.files' not in use
    if(length(unused)>0) {
      unused <- paste(unused, collapse=", ")
      warning(paste("unmatched pdb files:", unused))
    }

    ## Elements of 'ids' not in use
    if(length(nonmatch)>0) {
      nonmatch <- paste(nonmatch, collapse=", ")
      warning(paste("unmatched ids:", nonmatch))
    }
  }
  return(out)
}
