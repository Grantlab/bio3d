context("Testing correlation network analysis")

test_that("cna() and cnapath() work properly", {

  skip_on_cran()

  oops <- requireNamespace("igraph", quietly = TRUE)
  if (!oops) {
     skip('Need igraph installed to run this test')
  }

  ## Simple test with PDB ID 1HEL
  file <- system.file("examples/1hel.pdb",package="bio3d")
  capture.output(pdb <- read.pdb(file))
 
  ## Normal modes and cross-correlation matrix 
  capture.output(modes <- nma(pdb))
  capture.output(cij <- dccm(modes))

  ## network construction
  capture.output(cm <- cmap(pdb, dcut=4.5, scut=1))
  suppressWarnings(
    capture.output(net <- cna(cij, cm=cm, cutoff.cij=0))
  )
  
  expect_equal(net$communities$vcount, 129)
  expect_equal(igraph::ecount(net$network), 608)
  expect_equal(max(net$communities$membership), 5)
  expect_equal(igraph::ecount(net$community.network), 9)
 
  ## and path analysis
  capture.output(pa <- cnapath(net, from=3, to=53, k=20, collapse=TRUE, ncore=1))

  expect_equal(length(pa$path), 20)
  expect_equal(pa$path[[1]], c(3, 2, 39, 54, 53))
  expect_equal(pa$path[[20]], c(3,  2, 39, 40, 84, 53))
  expect_equal(pa$dist[1:6], c(3.427, 3.521, 3.821, 3.847, 3.868, 3.965), tolerance=1e-3)
  expect_equal(pa$epath[[10]], c(9,  13, 245, 253, 320))

  ndg0 <- matrix(c(0.65, 1, 0.35, 0.8, 0.45, 0.1, 0.3,  1, 0.5, 0.35), nrow=1,
                 dimnames = list(1, c(2, 3, 38, 39, 40, 41, 42, 53, 54, 55)) )
  capture.output(ndg <- summary(pa)$degeneracy)
  expect_equal(ndg, ndg0)

  skip_on_travis()

  capture.output(pa2 <- cnapath(net, from=3, to=53, k=20, collapse=TRUE, ncore=NULL))
  expect_equal(pa, pa2)
})

