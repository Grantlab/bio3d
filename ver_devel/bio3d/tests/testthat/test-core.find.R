context("Testing core.find function")


test_that("core.find() works properly", {
  skip_on_cran()

  invisible(capture.output(pdbfiles <- get.pdb(c("1bg2","2ncd","1i6i","1i5s"), URLonly=TRUE)))
  invisible(capture.output(pdbs <- pdbaln(pdbfiles)))

  ##-- Very rough 'core' finding: Just for test
  invisible(capture.output(core <- core.find(pdbs, stop.vol=100, ncore=1)))

  resnos <- c(7, 195,  45, 194, 193, 252, 192, 197, 198,
              44, 253, 150, 102, 220, 272, 151, 158, 157)
  expect_equal(length(core$resno), 294)
  expect_equal(resnos, as.numeric(core$resno[1:18]))
  

  ## Check multicore 
  invisible(capture.output(core.mc <- core.find(pdbs, stop.vol=100, ncore=NULL)))
  expect_identical(core, core.mc)
})
