context("Testing get.seq()")

test_that("get.seq() works properly", {
  skip_on_cran()
  skip_on_travis()

  outfile <- tempfile()
  
  # NR (use GI/RefSeq Accession)
  ids <- c("829581541", "NP_851365")
  capture.output( seqs <- get.seq(ids, outfile=outfile, db='nr') )
  expect_equal(length(seqs$id), 2)
  expect_identical(head(seqs$ali[1,]), c('M','E','L','E','N','I'))
  expect_identical(head(seqs$ali[2,]), c('M','G','A','G','A','S'))
  
  unlink(outfile)
  
  # uniprot/swissprot
  ids <- c("P11488", "P62873")
  capture.output(seqs <- get.seq(ids, outfile=outfile, db='uniprot'))
  expect_equal(length(seqs$id), 2)
  expect_identical(head(seqs$ali[1, ]), c('M','G','A','G','A','S'))
  expect_identical(head(seqs$ali[2, ]), c('M','S','E','L','D','Q'))
  
  unlink(outfile)
  
  # pdb
  ids <- c("1tag", "1tnd_B")
  capture.output(seqs <- get.seq(ids, outfile=outfile, db='pdb'))
  expect_equal(length(seqs$id), 2)
  expect_identical(head(seqs$ali[1,]), c('A','R','T','V','K','L'))
  expect_identical(tail(seqs$ali[2,]), c('K','D','C','G','L','F'))
  
  unlink(outfile)
  
  # one of the ids is wrong
  ids <- c("P11488", "P62873", "P999")
  expect_warning(capture.output(seqs <- get.seq(ids, outfile=outfile, db='uniprot')))
  expect_true(is.logical(seqs))
  expect_identical(as.vector(seqs), c(FALSE, FALSE, TRUE))
  
  unlink(outfile)
  
  # all ids are wrong (choose wrong database)
  ids <- c("1tag", "1tnd_B")
  expect_error( capture.output(seqs <- get.seq(ids, outfile=outfile, db='uniprot')) )
   
  unlink(outfile)
})
