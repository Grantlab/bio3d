
configuration <<- list(
  
  pdbdir = list(
    archive = FALSE,
    rawfiles = "raw_files",
    splitfiles = "split_files"
    ),
  
  db = list(
    use = FALSE,
    host = "localhost",
    dbname = "shiny",
    username = "shiny",
    password = ""
    ),
  
  hmmer = list(
    local = FALSE,
    exefile = system('which phmmer', intern=TRUE),
    pdbseq = "/path/to/pdb/seqs/pdb_seqres.txt"
    ),
  
  muscle = list(
    exefile = system('which muscle', intern=TRUE)
    ),
  
  user_data = list(
    path = "/tmp"
    )
  
  )
