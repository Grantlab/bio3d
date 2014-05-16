context("Testing aa2mass()")


test_that("Amino acid mass tests", {

  ## Simple test
  sequ   <- c("ALA", "LYS", "TPO")
  masses <- c(71.08, 129.184, 181.084)
  
  expect_that(aa2mass(sequ, addter=FALSE, mmtk=FALSE, mass.custom=NULL),
              equals(masses))

  ## With Terminal atoms added
  masses <- c(72.088, 129.184, 198.092)
  expect_that(aa2mass(sequ, addter=TRUE, mmtk=FALSE, mass.custom=NULL),
              equals(masses))
  
  ## With 'custom' residues
  sequ   <- c("MLY", "HMM", "UNK")
  masses <- c(156.228, 10, 20.001)
  expect_that(aa2mass(sequ, addter=FALSE, mmtk=FALSE,
                      mass.custom=list(HMM=10, UNK=20.001)),
              equals(masses))
  
  
  
  
}
          )
