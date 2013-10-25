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
}
)


test_that("struct.aln still works", {
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
   # transducin
   ids <- c("1tag_A", "1tnd_A", "3sn6_A", "1cul_C",
            "3v00_A", "1got_A", "1kjy_A", "4g5o_A",
            "2ode_A", "2ihb_A", "1as0_A", "1tad_A")
   tmp <- tempdir()
   raw.files <- get.pdb(substr(ids, 1, 4), tmp, verbose=FALSE)
   invisible(capture.output(files <- pdbsplit(raw.files, ids=ids, path=file.path(tmp, "split_chain"))))
   invisible(capture.output(pdbs <- pdbaln(files)))
   gaps.pos <- gap.inspect(pdbs$xyz)

   # check if ncore > 1 is really faster 
   time1 <- system.time(xyz1 <- fit.xyz(pdbs$xyz[1,], pdbs$xyz, gaps.pos$f.inds, gaps.pos$f.inds, ncore=1))
   time2 <- system.time(xyz2 <- fit.xyz(pdbs$xyz[1,], pdbs$xyz, gaps.pos$f.inds, gaps.pos$f.inds, ncore=2))
   time1 <- time1["elapsed"]
   time2 <- time2["elapsed"]

   expect_equivalent(xyz1, xyz2)

#   cat("Speed up by", round((time1-time2)/time2, 1)*100, "%", sep="")
#   if(getOption("cores") > 1)
#      expect_true(time1 > time2)
})

