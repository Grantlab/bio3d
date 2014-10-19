context("Testing pca()")

test_that("pca functions works", {

  attach(transducin)
  inds <- unlist(lapply(c("1TND_A", "1TAG", "1AS0", "1AS2"), grep, pdbs$id))
  pdbs <- pdbs.filter(pdbs, row.inds=inds)
  gaps <- gap.inspect(pdbs$xyz)

  ## Calc modes
  invisible(capture.output(pc <- pca(pdbs)))

  ## check dimensions
  expect_that(dim(pc$U), equals(c(939, 939)))
  expect_that(length(pc$L), equals(939))

  ## check values
  Lexpected <- c(1.964689e+02, 1.715903e+02, 7.091482e+01, 5.372565e-13, 4.100013e-13,
                 3.515269e-13, 2.677252e-13, 2.539279e-13, 2.497755e-13, 2.234660e-13)

  expect_that(head(pc$L, n=10), equals(Lexpected, tolerance=1e-6))

  AUexpected <- c(0.01342227, 0.02528470, 0.03448167, 0.23837239, 0.13487428, 0.13964307)
  expect_that(head(pc$au[1,], n=6), equals(AUexpected, tolerance=1e-6))

  Z1expected <- c(-3.555218e+00, -3.823427e+00, 1.220469e+01,
                  -1.327600e-14, -2.708901e-14, 1.247504e-14)

  expect_that(head(pc$z[1,], n=6), equals(Z1expected, tolerance=1e-6))

  Mexpected <- c(30.12193, 67.76449, 43.36594, 27.01919, 69.66411, 44.47434)
  expect_that(head(pc$mean, n=6), equals(Mexpected, tolerance=1e-6))

  detach(transducin)
})
