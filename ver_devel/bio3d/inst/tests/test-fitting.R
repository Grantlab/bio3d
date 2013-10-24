context("Testing fitting functions")


test_that("Fitting still works", {

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
  
  rmsd.mat <- matrix(c(0.000, 0.293, 0.286,
                       0.293, 0.000, 0.270,
                       0.286, 0.270, 0.000),
                     ncol=3, byrow=TRUE)
  expect_that(rmsd( xyz[, gaps$f.inds] ),
              equals(rmsd.mat, tolerance  = 1e-6))


  ## Test struct.aln
  invisible(capture.output(pdb.a <- read.pdb("1hel")))
  invisible(capture.output(pdb.b <- read.pdb("1dpx")))
  
  invisible(capture.output(aln <- struct.aln(pdb.a, pdb.b, write.pdbs=FALSE,
                                            cutoff=0.1, max.cycles=2)))
  rmsda <- c(0.293, 0.229, 0.200)

  
  expect_that(aln$rmsd, equals(rmsda, tolerance  = 1e-6))
  expect_that(length(aln$a.inds$atom), equals(112))
  expect_that(length(aln$b.inds$atom), equals(112))
  expect_that(length(aln$b.inds$xyz), equals(112*3))
  expect_that(length(aln$b.inds$xyz), equals(112*3))
  
}
)
