context("Testing nma.pdbs()")

test_that("eNMA works", {
  tmp <- tempdir()
  
  ids <- c("1CDK_A", "1CMK_E", "3DND_A")
  invisible(capture.output(raw.files <- get.pdb(ids, path = tmp, gzip=TRUE)))
  invisible(capture.output(files <- pdbsplit(raw.files, ids = ids, path = tmp)))
  invisible(capture.output(pdbs <- pdbaln(files, fit=TRUE)))

  ## Calc modes
  invisible(capture.output(modes <- nma.pdbs(pdbs, fit=TRUE, rm.gaps=TRUE, ncore=1)))

  ## check dimensions
  expect_that(dim(modes$U), equals(c(1008, 1002, 3)))
  expect_that(dim(modes$L), equals(c(3, 1002)))
  expect_that(dim(modes$fluctuations), equals(c(3, 336)))
  
  ## structure 1- mode1:
  U1 <- c(-0.01713449,  0.02317549,  0.02522696,
          -0.01565469,  0.02459163,  0.03433094)
  nowU1 <- head(modes$U.subspace[,1,1], n=6)
  expect_that(nowU1, equals(U1, tolerance=1e-6))

  ## structure 1- mode2:
  U2 <- c(0.022968348, -0.012008463, -0.008486024,
          0.014687499, -0.014436478, -0.005957756)
  nowU2 <- head(modes$U.subspace[,2,1], n=6)
  expect_that(nowU2, equals(U2, tolerance=1e-6))

  ## structure 2- mode3:
  U3 <- c(0.003026549,  0.031369659,  0.021926428,
          -0.004656332,  0.034181381, 0.036996015)
  nowU3 <- head(modes$U.subspace[,3,2], n=6)
  expect_that(nowU3, equals(U3, tolerance=1e-6))

  ## Fluctuations:
  f1 <- c(0.2799137, 0.3139167, 0.3421285, 0.2781802, 0.2188650, 0.2739127)
  f2 <- c(0.2600824, 0.3043772, 0.3363984, 0.2731422, 0.2240850, 0.2610135)
  f3 <- c(0.5796615, 0.5201151, 0.4703789, 0.3598020, 0.2312163, 0.2843633)
  expect_that(modes$fluctuations[1,1:6], equals(f1, tolerance=1e-6))
  expect_that(modes$fluctuations[2,1:6], equals(f2, tolerance=1e-6))
  expect_that(modes$fluctuations[3,1:6], equals(f3, tolerance=1e-6))
  
  ## Orthognal
  expect_that(as.numeric(modes$U.subspace[,1,1] %*% modes$U.subspace[,1,1]),
              equals(1))
  expect_that(as.numeric(modes$U.subspace[,1,1] %*% modes$U.subspace[,2,1]),
              equals(0, tolerance=1e-6))

  ## RMSIP
  rmsips <- c(1.0000, 0.8948,  0.9333,
              0.8948, 1.0000,  0.8899,
              0.9333, 0.8899,  1.0000)
              
  expect_that(c(modes$rmsip), equals(rmsips, tolerance=1e-6))

  ### Multicore (same arguments as above!)
  invisible(capture.output(mmc <- nma.pdbs(pdbs, fit=TRUE, rm.gaps=TRUE, ncore=4)))
  expect_that(mmc$fluctuations, equals(modes$fluctuations))
  expect_that(mmc$U.subspace, equals(modes$U.subspace))
  
  
})
          

  
