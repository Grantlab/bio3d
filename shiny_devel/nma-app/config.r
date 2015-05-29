## NOTE: dir must exist prior to starting server!
## check also permissions for the shiny user

configuration <<- list(

  pdbdir = list(
    archive = FALSE,
    #rawfiles = "/media/elephant/pdb/archive2",
    #splitfiles = "/media/elephant/pdb/archive2/split"
    
    rawfiles = "raw_files",
    splitfiles = "split_files"
    ),
  
  user_data = list(
    path = "/tmp/user_data"
    )
  
  )
