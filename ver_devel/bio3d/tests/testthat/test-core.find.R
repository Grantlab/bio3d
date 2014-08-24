context("Testing core.find function")


test_that("core.find() works properly", {

  invisible(capture.output(pdbfiles <- get.pdb(c("1bg2","2ncd","1i6i","1i5s"), URLonly=TRUE)))
  invisible(capture.output(pdbs <- pdbaln(pdbfiles)))

  ##-- Very rough 'core' finding: Just for test
  invisible(capture.output(core <- core.find(pdbs, stop.vol=100, ncore=1)))

  expect_equal(length(core$all.resno), 294)
  expect_equal(core$all.resno[1], "7")
  expect_equal(core$all.resno[2], "195")
  expect_equal(core$all.resno[3], "45")

  ## Check multicore 
  invisible(capture.output(core.mc <- core.find(pdbs, stop.vol=100, ncore=NULL)))
  expect_identical(core, core.mc)
})
