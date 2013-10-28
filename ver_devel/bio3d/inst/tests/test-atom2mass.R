context("Testing atom mass functions")


test_that("atom to mass tests", {

  ## Simple test
  atom.names <- c("CA", "O", "N", "OXT")
  masses <- c(12.01, 16.00, 14.01, 16.00)
  expect_that(atom2mass(atom.names), equals(masses))

  masses <- c(42.02, 16.00)
  expect_that(atom2mass(atom.names, grpby=c(1,1,1,2)), 
              equals(masses))
  
  ## Should end with error  
  atom.names <- c("CA", "O", "N", "OXT", "CL2", "PT1")
  expect_that(atom2mass(atom.names, rescue=FALSE), throws_error())
  expect_that(atom2mass(atom.names, rescue=TRUE), gives_warning())

  ## Try with a PDB object
  invisible(capture.output(pdb.small <- read.pdb("1etl")))
  expect_that(sum(atom2mass(pdb.small)), equals(1165.396))

  ## Try center of mass at the same go
  coma <- c(16.3714058,  5.0264399,  0.1739735)
  expect_that(com(pdb.small), equals(coma, tolerance=1e-6))
  
  
  ## Add custom masses
  elety.cust <- list("CL2"="Cl", "PT1"="Pt")
  mass.cust <- list("Cl"=35.45, "Pt"=195.08)
  atom.names <- c("CA", "O", "N", "OXT", "CL2", "PT1")
  masses <- c(12.01, 16.00, 14.01, 16.00, 35.45, 195.08)
  
  expect_that(atom2mass(atom.names, mass.custom=mass.cust, elety.custom=elety.cust),
              equals(masses))
  
  ## mass from formula
  form <- "C5 H6 N O3"
  masses <- c(60.050,  6.048, 14.010, 48.000)
  expect_that(formula2mass(form, sum.mass=FALSE),
              equals(masses))

  form <- "C5H6"
  expect_that(formula2mass(form), throws_error())

}
          )
