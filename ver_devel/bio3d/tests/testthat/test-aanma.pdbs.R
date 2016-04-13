context("Testing aanma.pdbs()")

test_that("aanma based eNMA works", {

  skip_on_cran()
  skip_on_travis()

  "mysign" <- function(a,b) {
    if(all(sign(a)==sign(b)))
      return(1)
    else
      return(-1)
  }
  
  pdbdir <- tempdir()
  invisible( capture.output(pdbfiles <- get.pdb(c("1TND_A", "1TAG", "1AS0", "1AS2"), path=pdbdir, split=TRUE)) )
  invisible( capture.output(pdbs <- pdbaln(pdbfiles)) )
  gaps <- gap.inspect(pdbs$xyz)

  ## Calc modes
  invisible(capture.output(modes <- aanma.pdbs(pdbs, ncore=1)))

  ## check dimensions
  expect_that(dim(modes$U), equals(c(939, 933, 4)))
  expect_that(dim(modes$L), equals(c(4, 933)))
  expect_that(dim(modes$fluctuations), equals(c(4, 313)))
  expect_that(dim(modes$rmsip), equals(c(4, 4)))
  
  ## structure 1- mode1:
  U1 <- c(0.04766153, -0.007993432, -0.06513697, 0.04313399, -0.008026397, -0.03755748)
  nowU1 <- head(modes$U.subspace[,1,1], n=6)
  expect_that(nowU1 * mysign(U1, nowU1), equals(U1, tolerance=1e-4))

  ## structure 1- mode2:
  U2 <- c(0.01151413, -0.01098311, -0.05809988, 0.006112359, -0.0164621, -0.04298916)
  nowU2 <- head(modes$U.subspace[,2,1], n=6)
  expect_that(nowU2 * mysign(U2, nowU2), equals(U2, tolerance=1e-4))

  ## structure 4- mode3:
  U3 <- c(-0.1188249, -0.05850153, 0.02082409, -0.09123394, -0.02270938, 0.00987617)
  nowU3 <- head(modes$U.subspace[,3,4], n=6)
  expect_that(nowU3 * mysign(U3, nowU3), equals(U3, tolerance=1e-4))

  ## structure 4-mode1 - tail:
  U1 <- c(0.009886232, -0.0006358995, -0.05189479, 0.009095127, 0.0005291993, -0.06936025)
  nowU1 <- tail(modes$U.subspace[,1,4], n=6)
  expect_that(nowU1 * mysign(U1, nowU1), equals(U1, tolerance=1e-4))


  ## Fluctuations:
  f1 <- c(0.4259852, 0.3173786, 0.2339003, 0.2126361, 0.1965271, 0.1865991)
  f2 <- c(0.5467513, 0.4003687, 0.2712613, 0.2266669, 0.2006127, 0.1855429)
  f4 <- c(0.8295543, 0.4094441, 0.281545, 0.2527504, 0.2091488, 0.1931585)

  expect_that(modes$fluctuations[1,1:6], equals(f1, tolerance=1e-3))
  expect_that(modes$fluctuations[2,1:6], equals(f2, tolerance=1e-3))
  expect_that(modes$fluctuations[4,1:6], equals(f4, tolerance=1e-3))
  
  ## Orthognal
  expect_that(as.numeric(modes$U.subspace[,1,1] %*% modes$U.subspace[,1,1]),
              equals(1, tolerance=1e-6))
  expect_that(as.numeric(modes$U.subspace[,1,1] %*% modes$U.subspace[,2,1]),
              equals(0, tolerance=1e-6))

  ## RMSIP
  rmsips <- c(1.0000, 0.8988234, 0.9627764, 0.915605)
  expect_that(as.vector(modes$rmsip[1,]), equals(rmsips, tolerance=1e-3))

  
  ## Multicore (same arguments as above!)
  invisible(capture.output(mmc <- aanma.pdbs(pdbs, ncore=NULL)))
  expect_that(mmc$fluctuations, equals(modes$fluctuations, tolerance=1e-6))
  expect_that(mmc$U.subspace, equals(modes$U.subspace, tolerance=1e-6))


  
  ## Calc modes with rm.gaps=FALSE
  invisible(capture.output(modes <- aanma.pdbs(pdbs, rm.gaps=FALSE, ncore=NULL)))

  ## structure 1-mode1 - tail:
  U1 <- c(-0.05024617, 0.009081938, 0.01764626, -0.09325497, 0.01754876, 0.01996309)
  nowU1 <- tail(modes$U.subspace[,1,1], n=6)
  expect_that(nowU1 * mysign(U1, nowU1), equals(U1, tolerance=1e-4))

  ## structure 2-mode1 - tail:
  U1 <- c(0.008412425, -0.005413275, -0.06655021, NA, NA, NA)
  nowU1 <- modes$U.subspace[940:945,1,2]
  U1[is.na(U1)] <- 0
  nowU1[is.na(nowU1)] <- 0
  expect_that(nowU1 * mysign(nowU1, U1), equals(U1, tolerance=1e-4))

  ## fluctuations
  na.expected <- c(3, 4, 1258, 1259, 1262, 1263, 1266, 1267, 1268, 1270, 1271, 1272, 1274, 1275, 1276,
                   1278, 1279, 1280, 1282, 1283, 1284, 1286, 1287, 1288, 1290, 1291, 1292)
                   
  expect_that(which(is.na(modes$fluctuations)), equals(na.expected))
  
 
  f1 <- c(0.9114199, 0.4010512, 0.2977837, 0.2202269)
  f4 <- c(0.4401748, 0.5743443, 0.6644086, rep(NA, 7))
          
  expect_that(modes$fluctuations[1,1:4], equals(f1, tolerance=1e-3))
  expect_that(tail(modes$fluctuations[4,], n=10), equals(f4, tolerance=1e-3))


  ## Calc modes with mass=FALSE and temp=NULL
  invisible(capture.output(modes <- aanma.pdbs(pdbs, mass=FALSE, temp=NULL, ncore=NULL)))

  ## structure 1- mode1:
  U1 <- c(0.04185427, -0.007185472, -0.05865138, 0.04704491, -0.009310229, -0.0432674)
  nowU1 <- head(modes$U.subspace[,1,1], n=6)
  expect_that(nowU1 * mysign(U1, nowU1), equals(U1, tolerance=1e-4))
  
  ## structure 1- mode2:
  U2 <- c(-0.006541777, 0.01006403, 0.04607164, -0.00293836, 0.01800689, 0.043626)
  nowU2 <- head(modes$U.subspace[,2,1], n=6)
  expect_that(nowU2 * mysign(U2, nowU2), equals(U2, tolerance=1e-4))

  ## structure 4- mode3:
  U3 <- c(-0.09327228, -0.04605804, 0.01065858, -0.08065877, -0.02255212, 0.005277435)
  nowU3 <- head(modes$U.subspace[,3,4], n=6)
  expect_that(nowU3 * mysign(U3, nowU3), equals(U3, tolerance=1e-4))

})
          

  
