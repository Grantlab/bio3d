context("Testing get.seq()")

test_that("get.seq() works properly", {
  skip_on_cran()
  skip_on_travis()

  outfile <- tempfile()
  
  # NR
  # Use mixed RefSeq, Uniprot, and PDB Accession.
  ids <- c("NP_851365", "P62873", "1L3R_E", "5P21_A")
  capture.output( seqs <- get.seq(ids, outfile=outfile, db='nr') )
  expect_equal(length(seqs$id), 4)
  expect_identical(head(seqs$ali[1,]), c('M','G','A','G','A','S'))
  expect_identical(head(seqs$ali[2,]), c('M','S','E','L','D','Q'))
  expect_identical(head(seqs$ali[3,]), c('G','N','A','A','A','A'))
  expect_identical(head(seqs$ali[4,]), c('M','T','E','Y','K','L'))
  unlink(outfile)
  
  Sys.sleep(5)
  
  # NR
  # Retrieved IDs may be different from the query; shoot warning
  # In the following, a GI number is converted to a PDB code.
  ids <- c("829581541")
  expect_warning( capture.output(seqs <- get.seq(ids, outfile=outfile, db='nr')) )
  expect_equal(length(seqs$id), 1)
  expect_identical(seqs$id, c("pdb|4TND|A"))
  expect_identical(head(seqs$ali[1,]), c('M','E','L','E','N','I'))
  unlink(outfile)
  
  Sys.sleep(5)
  
  # NR
  # Another example that retrieved ID changes is when query ID is from
  # a different database. 5X3X has both "A" and "a" chains in PDB.
  # But in 'nr', they are coded as "A" and "AA".
  # But the server behaves weirdly in this case:
  # - input '5X3X_A' or '5X3X_a' (bio3d calls toupper() internally) 
  #   gives either A or AA randomly.
  # - input '5X3X_AA' does not work anymore.
  # 
  # Comment out because of the random behavior
  #ids <- c("5X3X_A", "5X3X_a", "5X3X_AA")
  #expect_warning( seqs <- get.seq(ids, outfile=outfile, db='nr') )
  #expect_equal(length(seqs$id), 1)
  #expect_identical(seqs$id, c("pdb|5X3X|AA"))
  #expect_identical(head(seqs$ali[1,]), c('M','T','P','I','L','A'))
  #unlink(outfile)
  #Sys.sleep(5)

    
  # NR
  # Some wrong IDs.
  #         OK,       wrong chain ID,  typo
  ids <- c("4Q21_A",  "1TND_E",        "1L3_E")
  expect_warning( capture.output(seqs <- get.seq(ids, outfile=outfile, db='nr')) )
  expect_true(is.logical(seqs))
  expect_equivalent(seqs, c(FALSE, TRUE, TRUE))
  unlink(outfile)
  
  Sys.sleep(5)
  
  
  # EBI - uniprot/swissprot
  ids <- c("P11488", "P62873", "TCPA_YEAST")
  capture.output( seqs <- get.seq(ids, outfile=outfile, db='uniprot') )
  expect_equal(length(seqs$id), 3)
  expect_identical(head(seqs$ali[1, ]), c('M','G','A','G','A','S'))
  expect_identical(head(seqs$ali[2, ]), c('M','S','E','L','D','Q'))
  expect_identical(head(seqs$ali[3, ]), c('M','S','Q','L','F','N'))
  unlink(outfile)
  
  Sys.sleep(5)
  
  # EBI - pdb
  ids <- c("1TAG_A", "1TND_B", "3J3Q_0") 
  capture.output( seqs <- get.seq(ids, outfile=outfile, db='pdb') )
  expect_equal(length(seqs$id), 3)
  expect_identical(head(seqs$ali[1,]), c('A','R','T','V','K','L'))
  expect_identical(tail(seqs$ali[2,]), c('K','D','C','G','L','F'))
  expect_identical(head(seqs$ali[3,]), c('P','I','V','Q','N','L'))
  unlink(outfile)
  
  Sys.sleep(5)
  
  # EBI - pdb
  # Double-letter chain IDs are not supported.
  # So, for example, "3J3Q_10", "3J3Q_aA" etc. do not work.
  # 
  # ids <- c("3J3Q_10", "3J3Q_aA", "3J3Q_0")
  # expect_warning( seqs <- get.seq(ids, outfile=outfile, db='pdb') )
  # expect_true(is.logical(seqs))
  # expect_equivalent(seqs, c(TRUE, TRUE, FALSE))
  #
  # Smaller case chain IDs are not supported (because the internal 
  # case conversion in get.seq(); Does the server supports smaller cases?)
  #
  # Not all chains are available from the server:
  # ids <- c("3J3Q_A")
  # expect_error( seqs <- get.seq(ids, outfile=outfile, db='pdb') )
  #
  # In the following, the result will be constantly "A" (not random as NCBI)
  # ids <- c("5X3X_A")
  # seqs <- get.seq(ids, outfile=outfile, db='pdb') 
  #
  # Returns all capital cases only
  # ids <- c("5X3X")
  # seqs <- get.seq(ids, outfile=outfile, db='pdb') 
  
  
  # EBI - uniprot/swissprot
  # One of the ids is wrong.
  # Note: all sequences after the wrong id will not be retrieved.
  # The behavior for 'pdb' is different: only wrong ids are not retrieved.
  ids <- c("P11488", "P62873", "P999")
  expect_warning(capture.output(seqs <- get.seq(ids, outfile=outfile, db='uniprot')))
  expect_true(is.logical(seqs))
  expect_equivalent(seqs, c(FALSE, FALSE, TRUE))
  unlink(outfile)
  
  Sys.sleep(5)
  
  # EBI
  # All ids are wrong - choose wrong database
  ids <- c("1TAG_A", "1TND_B")
  expect_error( capture.output(seqs <- get.seq(ids, outfile=outfile, db='uniprot')) )
   
  unlink(outfile)
})
