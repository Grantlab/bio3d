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
  U7 <- c(-0.066793729, -0.053883268,  0.001891979,
          -0.036738304, -0.047565978,  0.001753576)
  nowU7 <- head(modes$U[,7])
  expect_that(nowU7 * mysign(U7, nowU7), equals(U7, tolerance=1e-4))

  ## Check second eigenvector
  U8 <- c(-0.0006820437, -0.0063015252, -0.0172583581,
          -0.0208773865, -0.0077078236, -0.0197741642)
  nowU8 <- head(modes$U[,8])
  expect_that(nowU8 * mysign(U8, nowU8), equals(U8, tolerance=1e-4))

  ## Check eigenvalues
  eival <- c(0.012147, 0.016651, 0.025490, 0.034212, 0.040073, 0.048434)
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
  U7 <- c(-0.059992473, -0.053860808,  0.002632462,
          -0.038492062, -0.060252657,  0.002212016)
  nowU7 <- head(modes$U[,7])
  expect_that(nowU7 * mysign(U7, nowU7), equals(U7, tolerance=1e-4))

  ## Check second eigenvector
  U8 <- c(-0.008837969, -0.003401558,  0.018160517, 
           0.018834639, -0.002279249,  0.022343009)
  nowU8 <- head(modes$U[,8])
  expect_that(nowU8 * mysign(U8, nowU8), equals(U8, tolerance=1e-4))

  ## Check eigenvalues
  eival <- c(0.005865, 0.007781, 0.010644, 0.015755, 0.019162, 0.020037)
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
  U7 <- c(0.067011206,  0.053922543, -0.001093240,
          0.036319837,  0.047227389, -0.001719102)
  nowU7 <- head(modes$U[,7])
  expect_that(nowU7 * mysign(U7, nowU7), equals(U7, tolerance=1e-4))

  ## Check second eigenvector
  U8 <- c(0.0007543764, 0.0061986220, 0.0164859315,
          0.0207461840, 0.0078183773, 0.0195714639)
  nowU8 <- head(modes$U[,8])
  expect_that(nowU8 * mysign(U8, nowU8), equals(U8, tolerance=1e-4))

  ## Check eigenvalues
  eival <- c(0.012450, 0.017159, 0.026374, 0.037100, 0.041601, 0.050559)
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
