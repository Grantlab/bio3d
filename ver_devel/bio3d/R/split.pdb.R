`split.pdb` <-
function(pdb.files, path="split_chain/") {
  ## split.pdb(pdb.files)
  toread <- file.exists(pdb.files)
  toread[substr(pdb.files, 1, 4) == "http"] <- TRUE
  if (all(!toread)) 
    stop("No corresponding PDB files found")

  if(any(!toread)) {
    warning(paste("Missing files:\n\t",
                  paste(pdb.files[!toread],
                        collapse="\n\t"), sep=""))
  
    pdb.files <- pdb.files[toread]
  }
  
  dir.create(path)
  for(i in 1:length(pdb.files)) {
    pdb <- read.pdb(pdb.files[i], het2atom=TRUE, maxlines=-1)
    chains <- unique(pdb$atom[,"chain"])
    if(length(chains) > 0) {
      for(j in 1:length(chains)) {
        if( !is.na(chains[j]) ) {
          new.pdb <- NULL
          
          sel <- paste("//",chains[j],"/////")
          new.pdb <- trim.pdb(pdb, atom.select(pdb,sel))
          
          new.name <- paste(substr(basename(pdb.files[i]),1,4),
                            "_", chains[j], ".pdb", sep="")
          new.name <- file.path( path, new.name)
          
          write.pdb(new.pdb, file=new.name)
        }
      }
    }
  }
}

