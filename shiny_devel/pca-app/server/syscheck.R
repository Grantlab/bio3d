
test_external <- function(exefile, arg="") {
  
  status <- system(paste(exefile, arg), ignore.stderr = TRUE, 
                   ignore.stdout = TRUE)
  
  if (!(status %in% c(0, 1)))
    return(FALSE)
  else 
    return(TRUE)
  
}

system_check <- function() {

  warns <- list()
  
  ## testing DSSP
  if(!test_external(configuration$dssp$exefile, "--version"))
    warns <- c(warns, "Launching DSSP failed")
  
  ## testing MUSCLE
  if(!test_external(configuration$muscle$exefile, "-version"))
    warns <- c(warns, "Launching MUSCLE failed")

  ## testing pmmer
  if(configuration$hmmer$local) {
    if(!test_external(configuration$hmmer$exefile, ""))
      warns <- c(warns, "Launching PHMMER failed")
    
    if(!file.exists(configuration$hmmer$pdbseq)) {
      warns <- c(warns, paste("File not found",
                              configuration$hmmer$pdbseq))
    }
  }

  ## testing database
  if(configuration$db$use) {
    con <- db_connect()

    if(is.logical(con)) {
      if(!con)
        warns <- c(warns, "Could not connect to MySQL database")

    }
    else {
      query <- "DESCRIBE pdb_annotation;"
      res <- dbGetQuery(con, query)
      expected <- "acc|chainId|chainLength|compound|db_id|experimentalTechnique|id|ligandId|ligandName|resolution|sequence|source|structureId"
      
      if(!paste(sort(res$Field), collapse="|") == expected)
        warns <- c(warns, "Missing fields in table 'pdb_annotation'")

      query <- "DESCRIBE pfam;"
      res <- dbGetQuery(con, query)
      
      expected <- "acc|chainId|eValue|id|pdbResNumEnd|pdbResNumStart|pfamAcc|pfamDesc|pfamName|structureId"
      
      if(!paste(sort(res$Field), collapse="|") == expected)
        warns <- c(warns, "Missing fields in table 'pfam'")
    }

    db_disconnect(con)
  }

  ## testing pdbdir
  if(!file.exists(configuration$pdbdir$rawfiles)) {
    warns <- c(warns, paste0("File not found: configuration$pdbdir$rawiles \n",
                             "   ('", configuration$pdbdir$rawfiles, "')"))
  }

  if(!file.exists(configuration$pdbdir$splitfiles)) {
    warns <- c(warns, paste0("File not found: configuration$pdbdir$splitfiles \n",
                             "   ('", configuration$pdbdir$splitfiles, "')"))
  }


  ## testing user_data dir
  if(!file.exists(configuration$user_data$path)) {
    warns <- c(warns, paste0("File not found: configuration$user_data$path \n",
                             "   ('", configuration$user_data$path, "')"))
  }

  if(length(warns) > 0) {
    message("\n\n###############################################\n")
    message("Warnings produced - check configuration file\n")
    lapply(warns, function(x) message(paste(" ", x)))
    message("\n#################################################")
    message("\n\n\n")
    return(FALSE)
  }
  else {
    return(TRUE)
  }
}




