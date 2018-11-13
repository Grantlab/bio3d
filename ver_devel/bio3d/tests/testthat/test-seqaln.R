context("Testing seqaln")


test_that("seqaln works", {
  skip_on_cran()
  skip_on_travis()

  ## seqaln with one sequence. should remove gaps
  seqs <- c("X", "-", "-", "A", "C", "A", "G", "K", "-")
  aln <- seqaln(seqs)
  expected <- c("X", "A", "C", "A", "G", "K")
  expect_identical(c(aln$ali), expected)
  

  ## align two sequences
  seqs <- seqbind(seqs,
                c("C", "A", "G", "G", "A", "G", "K"))
  aln <- seqaln(seqs)
  expected <- seqbind(c("-", "X", "A", "C", "A", "G", "K"),
                      c("C", "A", "G", "G", "A", "G", "K"))
  expect_identical(aln$ali, expected$ali)

    
  ## add a sequence to the (profile) alignment
  seq <- c("G", "A", "G", "K", "-")
  aln <- seqaln(seq, profile=aln)
  rownames(aln$ali) <- paste0("seq", 1:3)
  expected <- seqbind(c("-", "X", "A", "C", "A", "G", "K", "-"),
                      c("C", "A", "G", "G", "A", "G", "K", "-"),
                      c("-", "-", "-", "G", "A", "G", "K", "-"))
  expect_identical(aln$ali, expected$ali)

  ## test 'msa' option
  seqs <- get.seq(c("4q21_A", "1ftn_A"), outfile=tempfile())
  aln <- seqaln(seqs, outfile=tempfile())
  aln2 <- seqaln(seqs, outfile=tempfile(), exefile="msa")
  aln$call <- NULL; aln2$call <- NULL
  expect_identical(aln, aln2)

})

  
  
