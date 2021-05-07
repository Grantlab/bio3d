context("Testing relationship between read.cif and read.pdb")

test_that("read.cif() reads a normal pdb file", {

  skip_on_cran()
  skip_on_travis()
    
  ## Simple testing 
  file <- system.file("examples/1dpx.pdb",package="bio3d")
  invisible(capture.output(p1 <- read.pdb(file)))

  datdir <- tempdir()
  invisible(capture.output(get.pdb("1dpx", path=datdir, format="cif",
                                   overwrite = FALSE, verbose = FALSE)))
  
  suppressWarnings(
    invisible(capture.output(p2 <- read.cif(file.path(datdir, "1dpx.cif"))))
  )
  
  expect_is(p2$atom, "data.frame")
  expect_true(inherits(p2, "pdb"))
  expect_true(inherits(p2$xyz, "xyz"))
  expect_identical(p2$xyz, p1$xyz)

  expect_equal(nrow(p2$atom), nrow(p1$atom))
  expect_equal(sum(p2$calpha), sum(p1$calpha))
  expect_identical(p2$xyz, p1$xyz)
  expect_identical(p2$atom$type, p1$atom$type)
  ## offset for eleno between pdb and cif format
  #expect_identical(p2$atom$eleno, p1$atom$eleno)
  expect_identical(p2$atom$resid, p1$atom$resid)
  expect_identical(p2$atom$resno, p1$atom$resno)
  expect_identical(p2$atom$o, p1$atom$o)
  expect_identical(p2$atom$b, p1$atom$b)
  
  expect_equal(sum(p2$atom$resid=="HOH"), sum(p1$atom$resid=="HOH"))
  expect_equal(sum(p2$atom$resid=="CL"), sum(p1$atom$resid=="CL"))
  expect_that(sum(p2$xyz), equals(44657.12, tolerance=1e-6))

  expect_equal(sum(p2$atom$type=="ATOM"), sum(p1$atom$type=="ATOM"))
  expect_equal(sum(p2$atom$type=="HETATM"), sum(p1$atom$type=="HETATM"))

  
})


test_that("read.cif() on a multimodel object", {
    
  skip_on_cran()
  skip_on_travis()
  
  datdir <- tempdir()
  invisible(capture.output(get.pdb(c("1L2Y"), path=datdir,
                                   overwrite = FALSE, verbose = FALSE)))
  invisible(capture.output(get.pdb(c("1L2Y"), path=datdir, format="cif",
                                   overwrite = FALSE, verbose = FALSE)))
  # multi-model structure
  invisible(capture.output(p1 <- read.pdb(file.path(datdir, "1L2Y.pdb"), multi=TRUE)))
  suppressWarnings(
    invisible(capture.output(p2 <- read.cif(file.path(datdir, "1L2Y.cif"), multi=TRUE)))
  )
  
  expect_identical(dim(p1$xyz), dim(p2$xyz))
  expect_identical(p1$atom$x, p2$atom$x)
  expect_identical(p1$atom$y, p2$atom$y)
  expect_identical(p1$atom$z, p2$atom$z)
})
