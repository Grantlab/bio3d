context("Testing fitting functions")


test_that("Fitting still works", {
  skip_on_cran()

  ## Test fit.xyz / gap.inspect / pdbaln
  ##aln  <- read.fasta(system.file("examples/kif1a.fa",package="bio3d"))
  ##pdbs <- read.fasta.pdb(aln)
  
  invisible(capture.output(pdbs <- pdbaln(c("1hel", "1dpx", "1h87"))))
                                           
  gaps <- gap.inspect(pdbs$xyz)
   
  rmsd.mat <- matrix(c(0.000, 1.386, 1.245,
                       1.386, 0.000, 0.608,
                       1.245, 0.608, 0.000 ),
                     ncol=3, byrow=TRUE)

  ## Test rmsd()
  expect_that(rmsd( pdbs$xyz[, gaps$f.inds] ),
              equals(rmsd.mat, tolerance  = 1e-6))

  ## Test fit.xyz()
  xyz <- fit.xyz( fixed  = pdbs$xyz[1,],
                 mobile = pdbs$xyz,
                 fixed.inds  = gaps$f.inds,
                 mobile.inds = gaps$f.inds )

  x1.expected <- c(2.388, 10.533,  9.168,  2.524, 13.840,  7.282)
  x2.expected <- c(2.391083, 10.534633,  9.104952, 2.479335, 13.829487, 7.259277)
  x3.expected <- c(2.422049, 10.516869,  9.364834, 2.420416, 13.827526, 7.506205)
  expect_that(head(xyz[1,]), equals(x1.expected, tolerance  = 1e-6))
  expect_that(head(xyz[2,]), equals(x2.expected, tolerance  = 1e-6))
  expect_that(head(xyz[3,]), equals(x3.expected, tolerance  = 1e-6))

  rmsd.mat <- matrix(c(0.000, 0.293, 0.286,
                       0.293, 0.000, 0.270,
                       0.286, 0.270, 0.000),
                     ncol=3, byrow=TRUE)
  expect_that(rmsd( xyz[, gaps$f.inds] ),
              equals(rmsd.mat, tolerance  = 1e-6))

  xyz2 <- fit.xyz( fixed  = pdbs$xyz[1,],
                  mobile = pdbs,
                  fixed.inds  = gaps$f.inds,
                  mobile.inds = gaps$f.inds )
  expect_that(xyz2, equals(xyz))

}
)


test_that("struct.aln still works", {
  skip_on_cran()

  ## Test struct.aln
  invisible(capture.output(pdb.a <- read.pdb("1hel")))
  invisible(capture.output(pdb.b <- read.pdb("1dpx")))
  
  invisible(capture.output(aln <- struct.aln(pdb.a, pdb.b, write.pdbs=FALSE,
                                            cutoff=0.1, max.cycles=2, extra.args="-quiet")))
  rmsda <- c(0.293, 0.229, 0.200)

  
  expect_that(aln$rmsd, equals(rmsda, tolerance  = 1e-6))
  expect_that(length(aln$a.inds$atom), equals(112))
  expect_that(length(aln$b.inds$atom), equals(112))
  expect_that(length(aln$b.inds$xyz), equals(112*3))
  expect_that(length(aln$b.inds$xyz), equals(112*3))
}
)

# A little bit more tests...
test_that("fit.xyz() gets the same results as VMD", {
   skip_on_cran()

   invisible(capture.output(pdbs <- pdbaln(c("1tag", "1as0"))))
   inds <- gap.inspect(pdbs$xyz)$f.inds

   expect_error(fit.xyz("string", pdbs$xyz, inds, inds)) 
   expect_error(fit.xyz(pdbs$xyz, pdbs$xyz, inds, inds))
   expect_error(fit.xyz(pdbs$xyz[1,], pdbs$xyz, inds, 1:4))

   xyz <- fit.xyz(pdbs$xyz[1,], pdbs$xyz[2,], inds, inds)
   xyz <- xyz[!is.na(xyz)]
   xyz0 <- c(41.063, 15.667, 58.826, 43.041, 17.479, 56.113,
             44.826, 15.317, 53.571) # VMD results

   expect_equal(round(xyz[1:9], 3), xyz0)
})

test_that("fit.xyz() with ncore>1 works properly", {
   skip_on_cran()

   invisible(capture.output(pdbs <- pdbaln(c("1tag", "1as0", "1as2"))))
   inds <- gap.inspect(pdbs$xyz)$f.inds

   # check if ncore > 1 is really faster 
   time1 <- system.time(xyz1 <- fit.xyz(pdbs$xyz[1,], pdbs$xyz, inds, inds, ncore=1))
   time2 <- system.time(xyz2 <- fit.xyz(pdbs$xyz[1,], pdbs$xyz, inds, inds, ncore=NULL))
   time1 <- time1["elapsed"]
   time2 <- time2["elapsed"]

   expect_equivalent(xyz1, xyz2)

#   cat("Speed up by", round((time1-time2)/time2, 1)*100, "%", sep="")
#   if(getOption("cores") > 1)
#      expect_true(time1 > time2)
})

