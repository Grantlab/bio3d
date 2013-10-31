context("Testing deformation analysis")

test_that("still works", {
  
  invisible(capture.output(pdb.small <- read.pdb("1etl")))
  invisible(capture.output(modes <- nma(pdb.small)))
                           
  sums0 <- c(59.89283, 141.39431, 109.09525,
             122.52931, 172.63766, 317.01506)
  
  defe <- deformation.nma(modes)
  expect_that(defe$sums[1:6], equals(sums0, tolerance=1e-6))
  expect_that(defe$sums[1:6], equals(colSums(defe$ei[,1:6]), tolerance=1e-6))
   
})

test_that("fits with MMTK", {
  
  "calpha.mmtk" <- function(r, ...) {
    ## MMTK Units: kJ / mol / nm^2
    a <- 128; b <- 8.6 * 10^5; c <- 2.39 * 10^5;
    ifelse( r<4.0,
           b*(r/10) - c,
           a*(r/10)^(-6) )
  }

  ## Calc modes
  invisible(capture.output(pdb.small <- read.pdb("1etl")))
  invisible(capture.output(modes <- nma(pdb.small, pfc.fun=calpha.mmtk,
                                        addter=FALSE, mmtk=TRUE)))

  ## deformation energies of mode 7 (MMTK)
  def.mmtk <- c(1306.17014108, 524.571239022, 66.6665951865, 820.62710645,
                154.703500149, 754.482784094, 382.993752804, 173.118373857,
                287.880418213, 205.968139938, 466.277540766, 814.845931887)
  
  ## calc deformation energies
  defe <- deformation.nma(modes, mode.inds=seq(7,26), pfc.fun=calpha.mmtk)
  expect_that(defe$ei[,1], equals(def.mmtk, tolerance=1e-6))

  # mode 8
  def.mmtk <- c(1632.46638349, 4188.80097302, 1445.59255222, 1841.09326982)
  expect_that(head(defe$ei[,2], n=4), equals(def.mmtk, tolerance=1e-6))

  #mode 9
  def.mmtk <- c(1166.5762743, 917.219625799, 438.366075722, 1230.87278639)
  expect_that(head(defe$ei[,3], n=4), equals(def.mmtk, tolerance=1e-6))
   

})

