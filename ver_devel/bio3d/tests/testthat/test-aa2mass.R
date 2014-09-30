context("Testing aa2mass()")


test_that("Amino acid mass tests", {

  ## Simple test
  sequ   <- c("ALA", "LYS", "TPO")
  ## previous to Juliens atom mass lookup
  ##masses <- c(71.08, 129.184, 181.084)
  masses <- c(71.078, 129.180, 181.084)
  
  expect_that(aa2mass(sequ, addter=FALSE, mmtk=FALSE, mass.custom=NULL),
              equals(masses, tolerance=1e-6))

  ## With Terminal atoms added
  #masses <- c(72.088, 129.184, 198.092)
  masses <- c(72.08594, 129.18000, 198.09134)
  expect_that(aa2mass(sequ, addter=TRUE, mmtk=FALSE, mass.custom=NULL),
              equals(masses, tolerance=1e-6))
  
  ## With 'custom' residues
  sequ   <- c("MLY", "HMM", "UNK")
  masses <- c(156.225, 10.000, 20.001)
  expect_that(aa2mass(sequ, addter=FALSE, mmtk=FALSE,
                      mass.custom=list(HMM=10, UNK=20.001)),
              equals(masses, tolerance=1e-6))
   
}
          )
