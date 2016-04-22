context("Testing gnm()")


test_that("GNM", {

  "mysign" <- function(a,b) {
    if(all(sign(a)==sign(b)))
      return(1)
    else
      return(-1)
  }

  ## Simple test with PDB ID 1HEL
  file <- system.file("examples/1hel.pdb",package="bio3d")
  invisible(capture.output(pdb <- read.pdb(file)))
  
  ## Calculate modes with default arguments
  invisible(capture.output(modes <- gnm.pdb(pdb) ) )

  ## NOTE: all following reference results are from Prody
  ## Check first eigenvector 
  U2 <- c( -0.020,  0.008,  0.026,  0.052,  0.069,  0.086)
  nowU2 <- round(head(modes$U[, 2]), 3)
  expect_that(nowU2 * mysign(U2, nowU2), equals(U2, tolerance=1e-6))
  
  ## Check second eigenvector
  U3 <- c(0.050, 0.064, 0.084, 0.110, 0.105, 0.145)
  nowU3 <- round(head(modes$U[,3]), 3)
  expect_that(nowU3 * mysign(U3, nowU3), equals(U3, tolerance=1e-6))

  ## Check eigenvalues
  eival <- c(0.342,   0.804,   1.108,   1.277,   1.416,   1.617)
  nowEival <- round(modes$L[2:7], 3)
  expect_that(nowEival, equals(eival, tolerance=1e-6))

  ## Dimensions
  expect_that(dim(modes$U), equals(c(129, 129)))
  expect_that(length(modes$L), equals(129))
  expect_that(modes$natoms, equals(129))
  expect_that(modes$temp, equals(300))

  ## Orthognals
  expect_that(as.numeric(modes$U[,2] %*% modes$U[,2]),
              equals(1, tolerance=1e-6))
  expect_that(all(round(c(modes$U[,2] %*% modes$U[, 3:ncol(modes$U)]), 6)==0),
              equals(TRUE))
  
  expect_that(all(round(c(modes$L[1]), 6)==0), equals(TRUE))
  
  ## fluctuations (NOTE: Prody results are scaled here by the thermodynamic factor 3*k_B*T)
  flucts <- c( 1.379, 1.370, 1.108, 1.425, 1.044, 1.202, 1.315, 0.976, 0.855, 1.173 ) 
  nowFlucts <- round(modes$fluctuations[1:10], 3)
  expect_that(nowFlucts, equals(flucts, tolerance=1e-6))
  
  ## Covariance
  vcov <- c( 0.368,  0.293,  0.118,  0.032, -0.019,  0.037,  0.032, -0.040, -0.039, -0.020 )
  nowVcov <- round(cov.gnm(modes)[1, 2:11], 3)
  expect_that(nowVcov, equals(vcov, tolerance=1e-6))

 
}   )
