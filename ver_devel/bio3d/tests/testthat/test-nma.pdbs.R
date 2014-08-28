context("Testing nma.pdbs()")

test_that("eNMA works", {

  "mysign" <- function(a,b) {
    if(all(sign(a)==sign(b)))
      return(1)
    else
      return(-1)
  }

  tmp <- tempdir()
  
  ids <- c("1a70_A", "1czp_A", "1frd_A", "1fxi_A", "1iue_A", "1pfd_A")
  invisible(capture.output(raw.files <- get.pdb(ids, path = tmp, gzip=TRUE)))

  ##- previous to dataframe pdb format
  invisible(capture.output(files <- pdbsplit(raw.files, ids = ids, path = tmp)))
  invisible(capture.output(pdbs <- pdbaln(files, fit=TRUE)))

  ## Calc modes
  invisible(capture.output(modes <- nma.pdbs(pdbs, fit=TRUE, rm.gaps=TRUE, ncore=1)))

  ## check dimensions
  expect_that(dim(modes$U), equals(c(288, 282, 6)))
  expect_that(dim(modes$L), equals(c(6, 282)))
  expect_that(dim(modes$fluctuations), equals(c(6, 96)))
  
  ## structure 1- mode1:
  U1 <- c(-0.06121500,  0.10505209, -0.05435506, -0.04783683, 0.06649716, -0.03675693)
  nowU1 <- head(modes$U.subspace[,1,1], n=6)
  expect_that(nowU1 * mysign(U1, nowU1), equals(U1, tolerance=1e-6))

  ## structure 1- mode2:
  ##- previous to fixing mass calculation in nma.pdbs
  U2 <- c(0.038119084, -0.024193962, -0.025397837, 0.032008760, -0.005760642, -0.026255866)
  nowU2 <- head(modes$U.subspace[,2,1], n=6)
  expect_that(nowU2 * mysign(U2, nowU2), equals(U2, tolerance=1e-6))

  ## structure 4- mode3:
  ##- previous to fixing mass calculation in nma.pdbs
  U3 <- c(0.029950724, -0.004793673, -0.075769300, 0.039328324, 0.003828264, -0.059507805)
  nowU3 <- head(modes$U.subspace[,3,4], n=6)
  expect_that(nowU3 * mysign(U3, nowU3), equals(U3, tolerance=1e-6))

  ## Fluctuations:
  ##- previous to fixing mass calculation in nma.pdbs
  f1 <- c(0.44800568, 0.13568048, 0.11483728, 0.09895198, 0.10277404, 0.09565510)
  f2 <- c(0.38784015, 0.13387508, 0.10242378, 0.09160573, 0.09356470, 0.09255706)
  f6 <- c(0.28144615, 0.11132125, 0.09728422, 0.08709358, 0.08404063, 0.08160931)
  expect_that(modes$fluctuations[1,1:6], equals(f1, tolerance=1e-6))
  expect_that(modes$fluctuations[2,1:6], equals(f2, tolerance=1e-6))
  expect_that(modes$fluctuations[6,1:6], equals(f6, tolerance=1e-6))
  
  ## Orthognal
  expect_that(as.numeric(modes$U.subspace[,1,1] %*% modes$U.subspace[,1,1]),
              equals(1))
  expect_that(as.numeric(modes$U.subspace[,1,1] %*% modes$U.subspace[,2,1]),
              equals(0, tolerance=1e-6))

  ## RMSIP
  rmsips <- c(1.0000, 0.8355, 0.8788, 0.8348, 0.8977, 0.8598)
  expect_that(as.vector(modes$rmsip[1,]), equals(rmsips, tolerance=1e-6))


  
  ## Multicore (same arguments as above!)
  invisible(capture.output(mmc <- nma.pdbs(pdbs, fit=TRUE, rm.gaps=TRUE, ncore=3)))
  expect_that(mmc$fluctuations, equals(modes$fluctuations))
  expect_that(mmc$U.subspace, equals(modes$U.subspace))


  
  ## Calc modes with rm.gaps=FALSE
  invisible(capture.output(modes <- nma.pdbs(pdbs, fit=TRUE, rm.gaps=FALSE, ncore=3)))
  na.expected <- c(43, 46, 47, 48, 49, 52, 53, 54, 590,
                   591, 592, 594, 595, 596, 597, 598, 600)
  expect_that(which(is.na(modes$fluctuations)), equals(na.expected))
  
 
  f1 <- c(0.449953, 0.13691985, 0.11566978, 0.09961, 0.10306614,
          0.09598316, 0.10639936, NA, NA, 0.11016618)
  f2 <- c(0.38827622, 0.13437088, 0.10290902, 0.09219309, 0.0937538,
          0.09217573, 0.09791753, 0.10069342, 0.13709235, 0.31764473)
          
  f6 <- c(0.06942637, 0.08623632, 0.1085492, 0.08020983, 0.10722807,
          0.19175073, 0.1840169, 0.20978268, 0.33441566, NA, NA)
          
  expect_that(modes$fluctuations[1,1:10], equals(f1, tolerance=1e-6))
  expect_that(modes$fluctuations[2,1:10], equals(f2, tolerance=1e-6))
  expect_that(modes$fluctuations[6,90:100], equals(f6, tolerance=1e-6))
 
  
})
          

  
