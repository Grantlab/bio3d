context("Testing RMSD function")

test_that("rmsd() gets the same results as VMD", {
   invisible(capture.output(pdbs <- pdbaln(c("1tag", "1as0"))))
   inds <- gap.inspect(pdbs$xyz)$f.inds

   rmsd <- rmsd(pdbs$xyz[1, ], pdbs$xyz[2, ], inds, inds, fit=TRUE, ncore=1)
   rmsd0 <- 1.659 # VMD results

   expect_equal(round(rmsd, 3), rmsd0)
})

test_that("rmsd() with ncore>1 works properly", {
   invisible(capture.output(pdbs <- pdbaln(c("1tag", "1as0", "1as2"))))
   inds <- gap.inspect(pdbs$xyz)$f.inds

   # check if ncore > 1 is really faster 
   time1 <- system.time(rmsd1 <- rmsd(pdbs$xyz, a.inds=inds, ncore=1))
   time2 <- system.time(rmsd2 <- rmsd(pdbs$xyz, a.inds=inds, ncore=NULL))
   time1 <- time1["elapsed"]
   time2 <- time2["elapsed"]

   expect_equivalent(rmsd1, rmsd2)

#   cat("Speed up by", round((time1-time2)/time2, 1)*100, "%", sep="")
#   if(getOption("cores") > 1)
#      expect_true(time1 > time2)
})

