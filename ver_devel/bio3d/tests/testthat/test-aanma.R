context("Testing aanma()")

test_that("aanma based all-atom NMA works", {

  skip_on_cran()
#  skip_on_travis()

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
  invisible(capture.output(modes <- aanma(pdb, pfc.fun=load.enmff('aaenm2'),
                                        mass=TRUE, temp=300.0)))
 
  ## Check first eigenvector
  U7 <- c(0.06866314, 0.05545397, -0.00219604,
          0.03743162, 0.04853997, -0.001929003)
  nowU7 <- head(modes$U[,7])
  expect_that(nowU7 * mysign(U7, nowU7), equals(U7, tolerance=1e-4))

  ## Check second eigenvector
  U8 <- c(-0.0004398378, -0.006375451, -0.01791086, 
          -0.02096124,   -0.007711646, -0.02020036 )
  nowU8 <- head(modes$U[,8])
  expect_that(nowU8 * mysign(U8, nowU8), equals(U8, tolerance=1e-4))

  ## Check eigenvalues
  eival <- c(0.011362, 0.015574, 0.023757, 0.03164, 0.03751, 0.045251)
  nowEival <- modes$L[7:12]
  expect_that(nowEival, equals(eival, tolerance=1e-4))

  ## Dimensions
  expect_that(dim(modes$U), equals(c(387, 387)))
  expect_that(dim(modes$modes), equals(c(387, 387)))
  expect_that(length(modes$L), equals(387))
  expect_that(length(modes$frequencies), equals(387))
  expect_that(length(modes$mass), equals(129))
  expect_that(modes$natoms, equals(129))
  expect_that(modes$temp, equals(300))

  ## Orthognals
  tmpU <- crossprod(modes$U, modes$U)
  expect_true(all(round(diag(tmpU), 6)==1))
  expect_true(all(round(tmpU[upper.tri(tmpU)], 6)==0))

  expect_that(all(round(c(modes$L[1:6]), 6)==0), equals(TRUE))


  ###################################################################
  #
  # Test with reduced = TRUE 
  #
  ###################################################################

  ## Calculate modes
  invisible(capture.output(modes <- aanma(pdb, pfc.fun=load.enmff('aaenm2'),
                                        mass=TRUE, temp=300.0, reduced=TRUE)))
 
  ## Check first eigenvector
  U7 <- c(0.07243077, 0.06455693, -0.004789346, 
          0.03705018, 0.05790539, -0.00305091)
  nowU7 <- head(modes$U[,7])
  expect_that(nowU7 * mysign(U7, nowU7), equals(U7, tolerance=1e-4))

  ## Check second eigenvector
  U8 <- c(-0.01133821, -0.003806995, 0.0262246, 
           0.01619017, -0.002352164, 0.02449826)
  nowU8 <- head(modes$U[,8])
  expect_that(nowU8 * mysign(U8, nowU8), equals(U8, tolerance=1e-4))

  ## Check eigenvalues
  eival <- c(0.003232, 0.004269, 0.005667, 0.008636, 0.009515, 0.010767)
  nowEival <- modes$L[7:12]
  expect_that(nowEival, equals(eival, tolerance=1e-4))

  ## Orthognals
  tmpU <- crossprod(modes$U, modes$U)
  expect_true(all(round(diag(tmpU), 6)==1))
  expect_true(all(round(tmpU[upper.tri(tmpU)], 6)==0))
  
  expect_that(all(round(c(modes$L[1:6]), 6)==0), equals(TRUE))

 
  ###################################################################
  #
  # Test with rtb = TRUE 
  #
  ###################################################################

  ## Calculate modes
  invisible(capture.output(modes <- aanma(pdb, pfc.fun=load.enmff('aaenm2'),
                                        mass=TRUE, temp=300.0, rtb=TRUE)))
 
  ## Check first eigenvector
  U7 <- c(0.06890201, 0.05549496, -0.001374399, 
          0.03701697, 0.04819402, -0.001895817)
  nowU7 <- head(modes$U[,7])
  expect_that(nowU7 * mysign(U7, nowU7), equals(U7, tolerance=1e-4))

  ## Check second eigenvector
  U8 <- c(0.0005610636, 0.006274672, 0.0171175, 
          0.0208543, 0.007827309, 0.01999255)
  nowU8 <- head(modes$U[,8])
  expect_that(nowU8 * mysign(U8, nowU8), equals(U8, tolerance=1e-4))

  ## Check eigenvalues
  eival <- c(0.011646, 0.016051, 0.024595, 0.034361, 0.038952, 0.047275)
  nowEival <- modes$L[7:12]
  expect_that(nowEival, equals(eival, tolerance=1e-4))

  ## Orthognals
  tmpU <- crossprod(modes$U, modes$U)
  expect_true(all(round(diag(tmpU), 6)==1))
  expect_true(all(round(tmpU[upper.tri(tmpU)], 6)==0))
  
  expect_that(all(round(c(modes$L[1:6]), 6)==0), equals(TRUE))

  
  ###################################################################
  #
  # Test with outmodes = 'noh' 
  #
  ###################################################################

  ## Calculate modes
  invisible(capture.output(modes <- aanma(pdb, pfc.fun=load.enmff('aaenm2'),
                                        mass=TRUE, temp=300.0, outmodes='noh')))
 
  ## Check first eigenvector
  U7 <- c(0.0218632, 0.02019693, -0.002396826, 
          0.01873891, 0.01699585, -0.001568454)
  nowU7 <- head(modes$U[,7])
  expect_that(nowU7 * mysign(U7, nowU7), equals(U7, tolerance=1e-4))

  ## Check second eigenvector
  U8 <- c(0.009318529, 0.001879826, -0.005183737, 
          0.007078334, 0.0008036141, -0.003830362)
  nowU8 <- head(modes$U[,8])
  expect_that(nowU8 * mysign(U8, nowU8), equals(U8, tolerance=1e-4))

  ## Check eigenvalues
  eival <- c(0.011278, 0.013931, 0.01582, 0.023749, 0.030427, 0.039051)
  nowEival <- modes$L[7:12]
  expect_that(nowEival, equals(eival, tolerance=1e-4))

  ## Dimensions
  expect_that(dim(modes$U), equals(c(3003, 3003)))
  expect_that(dim(modes$modes), equals(c(3003, 3003)))
  expect_that(length(modes$L), equals(3003))
  expect_that(length(modes$frequencies), equals(3003))
  expect_that(length(modes$mass), equals(1001))
  expect_that(modes$natoms, equals(1001)) 
  expect_that(modes$temp, equals(300))

  ## Orthognals
  tmpU <- crossprod(modes$U[, 1:300], modes$U[, 1:300])
  expect_true(all(round(diag(tmpU), 6)==1))
  expect_true(all(round(tmpU[upper.tri(tmpU)], 6)==0))
  
  expect_that(all(round(c(modes$L[1:6]), 6)==0), equals(TRUE))

    
  ###################################################################
  #
  # Test with outmodes = 'noh' and rtb = TRUE
  #
  ###################################################################

  ## Calculate modes
  invisible(capture.output(modes <- aanma(pdb, pfc.fun=load.enmff('aaenm2'),
                               mass=TRUE, temp=300.0, outmodes='noh', rtb=TRUE)))
 
  ## Check first eigenvector
  U7 <- c(-0.0229911, -0.02044615, 0.002491489, 
          -0.01944523, -0.01711079, 0.001397449)
  nowU7 <- head(modes$U[,7])
  expect_that(nowU7 * mysign(U7, nowU7), equals(U7, tolerance=1e-4))

  ## Check second eigenvector
  U8 <- c(-0.005763786, -8.247862e-06, 0.00534345, 
          -0.003811432, 0.000763913, 0.003948311)
  nowU8 <- head(modes$U[,8])
  expect_that(nowU8 * mysign(U8, nowU8), equals(U8, tolerance=1e-4))

  ## Check eigenvalues
  eival <- c(0.011579, 0.014595, 0.016808, 0.024449, 0.032666, 0.040224)
  nowEival <- modes$L[7:12]
  expect_that(nowEival, equals(eival, tolerance=1e-4))

  ## Orthognals
  tmpU <- crossprod(modes$U[, 1:300], modes$U[, 1:300])
  expect_true(all(round(diag(tmpU), 6)==1))
  expect_true(all(round(tmpU[upper.tri(tmpU)], 6)==0))
  
  expect_that(all(round(c(modes$L[1:6]), 6)==0), equals(TRUE))

})
