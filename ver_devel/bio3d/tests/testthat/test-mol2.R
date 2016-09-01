context("Testing basic PDB structure operation")

test_that("read.mol2() reads a mol2 from zinc", {

  ## Simple test with aspirin
  file <- system.file("examples/aspirin.mol2",package="bio3d")
  invisible(capture.output(mol <- read.mol2(file)))

  expect_is(mol$atom, "data.frame")
  expect_true(inherits(mol, "mol2"))
  expect_true(inherits(mol$xyz, "xyz"))
  expect_equal(nrow(mol$atom), 20)
  expect_equal(nrow(mol$bond), 20)
  expect_equal(mol$info[1], 20)
  expect_equal(mol$info[2], 20)
  expect_equal(sum(mol$atom$elety=="H"), 7)

  elena <- c("C1", "C2", "O1", "O2", "C3", "C4", "C5",
             "C6", "C7", "C8", "C9", "O3", "O4", "H1", "H2",
             "H3", "H4", "H5", "H6", "H7")
  expect_equal(mol$atom$elena, elena)
      
  x <- c(-1.4238, -1.3441, -1.4532, -1.1519, -0.9822)
  y <- c(-2.5790, -3.9491, -4.7933, -4.2739, -2.8882)
  z <- c(0.6434,  0.0976,  0.0938,  0.0032, -0.0844)
  expect_equal(mol$atom$x[1:5], x)
  expect_equal(mol$atom$y[6:10], y)
  expect_equal(mol$atom$z[16:20], z)
  
})


test_that("read.mol2() reads and stores data properly", {
  skip_on_cran()
  skip_on_travis()
  
  file <- system.file("examples/aspirin.mol2",package="bio3d")
  invisible(capture.output(mol <- read.mol2(file)))

  f <- tempfile()
  write.mol2(mol, file=f)
  mol2 <- read.mol2(f)
  expect_equal(mol, mol2)
  
})


test_that("basic atom select and trim of mol2", {
  skip_on_cran()
  skip_on_travis()
  
  file <- system.file("examples/aspirin.mol2",package="bio3d")
  invisible(capture.output(mol <- read.mol2(file)))

  sele <- atom.select(mol, "noh")
  expect_equal(length(sele$atom), 13)

  sele <- atom.select(mol, elety="H")
  expect_equal(length(sele$atom), 7)

  sele <- atom.select(mol, elena="C1")
  expect_equal(length(sele$atom), 1)

  sele <- atom.select(mol, resno=1)
  expect_equal(length(sele$atom), 20)

  sele <- atom.select(mol, "noh")
  mol2 <- trim(mol, sele)

  expect_equal(nrow(mol2$atom), 13)
  expect_equal(nrow(mol2$bond), 13)
  expect_equal(length(mol2$xyz), 39)
  
  xyz <- c(-1.4238,  1.4221,  1.2577, -1.3441, -0.0813)
  expect_equal(mol2$xyz[1:5], xyz)
  
})



test_that("converting mol2 to pdb works", {
  skip_on_cran()
  skip_on_travis()
  
  file <- system.file("examples/aspirin.mol2",package="bio3d")
  invisible(capture.output(mol <- read.mol2(file)))

  pdb <- as.pdb(mol)
  expect_equal(nrow(pdb$atom), nrow(mol$atom))
  expect_equal(pdb$xyz, mol$xyz)

  expect_equal(mol$atom$elena, pdb$atom$elety)
  expect_equal(mol$atom$x, pdb$atom$x)
  expect_equal(mol$atom$charge, pdb$atom$charge)

})
