context("Testing pca()")

test_that("pca functions works", {
  tmp <- tempdir()
  
  ids <- c("1a70_A", "1czp_A", "1frd_A", "1fxi_A", "1iue_A", "1pfd_A")
  invisible(capture.output(raw.files <- get.pdb(ids, path = tmp, gzip=TRUE)))

  ##- previous to dataframe pdb format
  invisible(capture.output(files <- pdbsplit(raw.files, ids = ids, path = tmp)))
  invisible(capture.output(pdbs <- pdbaln(files, fit=TRUE)))
  
  ## Calc modes
  invisible(capture.output(pc <- pca(pdbs)))

  ## check dimensions
  expect_that(dim(pc$U), equals(c(288, 288)))
  expect_that(length(pc$L), equals(288))

  ## check values
  Lexpected <- c(5.458694e+01, 2.928604e+01, 1.249129e+01, 8.870719e+00,
                 4.535070e+00, 4.469946e-14, 4.255696e-14, 3.548377e-14,
                 3.067546e-14, 2.608522e-14)

  expect_that(head(pc$L, n=10), equals(Lexpected, tolerance=1e-6))

  AUexpected <- c(0.07777464, 0.06359141, 0.10253198, 0.05635136, 0.03189214, 0.12675553)
  expect_that(head(pc$au[1,], n=6), equals(AUexpected, tolerance=1e-6))

  Z1expected <- c(1.719149e+00, -2.258608e+00, -4.247174e+00,
                  2.319438e+00,  2.925958e+00, 6.578071e-15)

  expect_that(head(pc$z[1,], n=6), equals(Z1expected, tolerance=1e-6))

  Mexpected <- c(2.607529, -4.983709,  5.411727,  6.313369, -4.391827,  5.265203)
  expect_that(head(pc$mean, n=6), equals(Mexpected, tolerance=1e-6))
})
