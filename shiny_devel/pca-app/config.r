## wget ftp://ftp.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt.gz

configuration <<- list(

  pdbdir = list(
    rawfiles = "/path/to/pdb/raw_files",
    splitfiles = "/path/to/pdb/split_files"
    ),
  
  db = list(
    use = FALSE,
    host = "localhost",
    dbname = "test",
    username = "root",
    password = ""
    ),

  hmmer = list(
    local = FALSE,
    exefile = "phmmer",
    pdbseq = "/path/to/pdb/seqs/pdb_seqres.txt"
    ),

  muscle = list(
    exefile = "/usr/local/bin/muscle"
    )

  )
