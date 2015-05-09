## wget ftp://ftp.wwpdb.org/pub/pdb/derived <- data/pdb <- seqres.txt.gz

configuration <<- list(

                         pdbdir = list(
                                           rawfiles = "raw_files",
                                               splitfiles = "split_files"
                                               ),

                         db = list(
                                       use = FALSE,
                                           host = "localhost",
                                           dbname = "shiny_pdb-anno",
                                               username = "root",
                                               password = ""
                                                   ),

                         hmmer = list(
                                          local = FALSE,
                                              exefile = "phmmer",
                                              pdbseq = "/path/to/pdb/seqs/pdb_seqres.txt"
                                                  ),

                         muscle = list(
                                           exefile = system('which muscle', intern=TRUE)
                                               ),

                         user <- data = list(
                                                 path = "/tmp"
                                                     )

                           )
