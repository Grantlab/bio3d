context("Testing nma.pdbs()")

test_that("eNMA works", {
  tmp <- tempdir()
  
  ids <- c("1CDK_A", "1CMK_E", "3DND_A")
  invisible(capture.output(raw.files <- get.pdb(ids, path = tmp, gzip=TRUE)))
  invisible(capture.output(files <- pdbsplit(raw.files, ids = ids, path = tmp, 
                    het2atom=TRUE)))
  invisible(capture.output(pdbs <- pdbaln(files, fit=TRUE)))

  ## Calc modes
  invisible(capture.output(modes <- nma.pdbs(pdbs, fit=TRUE, rm.gaps=TRUE, ncore=1)))

  ## check dimensions
  expect_that(dim(modes$U), equals(c(1008, 20, 3)))
  expect_that(dim(modes$L), equals(c(3, 20)))
  expect_that(dim(modes$fluctuations), equals(c(3, 336)))
  
  ## structure 1- mode1:
  U1 <- c(0.01786246, -0.02389142, -0.02582718,
          0.01604795, -0.02540172, -0.03498606)
  nowU1 <- head(modes$U.subspace[,1,1], n=6)
  expect_that(nowU1, equals(U1, tolerance=1e-6))

  ## structure 1- mode2:
  U2 <- c(-0.022345000,  0.011023258,  0.007367781,
          -0.014114932,  0.013401786,  0.004408638)
  nowU2 <- head(modes$U.subspace[,2,1], n=6)
  expect_that(nowU2, equals(U2, tolerance=1e-6))

  ## structure 2- mode3:
  U3 <- c(-0.003280028, -0.031958422, -0.022230036,
          0.004565429,  -0.034818546, -0.037598946)
  nowU3 <- head(modes$U.subspace[,3,2], n=6)
  expect_that(nowU3, equals(U3, tolerance=1e-6))

  ## Fluctuations:
  f1 <- c(0.2798631, 0.3138516, 0.3420618, 0.2781450, 0.2188392, 0.2738726)
  f2 <- c(0.2598829, 0.3041458, 0.3361500, 0.2729427, 0.2239635, 0.2608458)
  f3 <- c(0.5795548, 0.5200072, 0.4702284, 0.3596465, 0.2311072, 0.2841777)
  expect_that(modes$fluctuations[1,1:6], equals(f1, tolerance=1e-6))
  expect_that(modes$fluctuations[2,1:6], equals(f2, tolerance=1e-6))
  expect_that(modes$fluctuations[3,1:6], equals(f3, tolerance=1e-6))
  
  ## Orthognal
  expect_that(as.numeric(modes$U.subspace[,1,1] %*% modes$U.subspace[,1,1]),
              equals(1))
  expect_that(as.numeric(modes$U.subspace[,1,1] %*% modes$U.subspace[,2,1]),
              equals(0, tolerance=1e-6))

  ## RMSIP
  rmsips <- c(1.0000, 0.8971, 0.9353,
              0.8971, 1.0000, 0.8870,
              0.9353, 0.8870, 1.0000)
  expect_that(c(modes$rmsip), equals(rmsips, tolerance=1e-6))

  ### Multicore (same arguments as above!)
  invisible(capture.output(mmc <- nma.pdbs(pdbs, fit=TRUE, rm.gaps=TRUE, ncore=4)))
  expect_that(mmc$fluctuations, equals(modes$fluctuations))
  expect_that(mmc$U.subspace, equals(modes$U.subspace))
  
  
})
          

  
