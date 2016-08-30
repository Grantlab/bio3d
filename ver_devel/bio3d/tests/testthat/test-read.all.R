context("Test reading aligned PDB structures")

test_that("read.all() reads all/select PDB atoms properly", {

  ## Internet access required - test skipped
  skip_on_cran()
  skip_on_travis()

  ## Test with 4 G-alpha PDB structures  
  pdbdir <- tempdir()
  invisible( capture.output(pdbfiles <- get.pdb(c("1TND_A", "1TAG", "1AS0", "1AS2"), path=pdbdir, split=TRUE)) )
  invisible( capture.output(aln <- pdbaln(pdbfiles)) )

  ## read with default options
  invisible( capture.output(pdbs <- read.all(aln, ncore=1)) )

  expect_is(pdbs, "pdbs")
  expect_true(inherits(pdbs, "fasta"))
  expect_is(pdbs$xyz, "xyz")
  expect_is(pdbs$all, "xyz")
  expect_equal(length(pdbs$all.hetatm), 4)
  expect_is(pdbs$all.hetatm[[1]], 'pdb')
  expect_equal(nrow(pdbs$all.hetatm[[1]]$atom), 35)

  xyz0 <- rbind(c(30.074,  64.434, 42.910, 30.390, 68.173, 43.382),
                c(38.552,  16.715, 60.576, 41.952, 15.364, 59.564),
                c(NA,      NA,     NA,     8.994, -26.463, 7.153),
                c(NA,      NA,     NA,     55.761, 11.671, 39.809))
  expect_equivalent(xyz0, pdbs$xyz[, 1:6])

  all0 <- rbind(c(30.848,  63.431,  43.718, 30.074,  64.434,  42.910),
                c(38.238,  18.018,  61.225, 38.552,  16.715,  60.576), 
                c(NA,      NA,     NA,     NA,   NA,   NA),
                c(NA,      NA,     NA,     NA,   NA,   NA))
  expect_equivalent(all0, pdbs$all[, 1:6])

  resno0 <- rbind(c(27, 28, 29, 30, 31, 32),
                  c(27, 28, 29, 30, 31, 32),
                  c(NA, 32, 33, 34, 35, 36),
                  c(NA, 32, 33, 34, 35, 36))
  expect_equivalent(resno0, pdbs$resno[, 1:6])

  b0 <- rbind(c(48.05, 40.48, 38.06, 29.55, 25.63, 23.82), 
              c(40.32, 31.04, 23.50, 15.42, 13.10, 13.73),
              c(NA,    51.59, 45.39, 31.58, 29.26, 22.08),
              c(NA,    88.56, 65.25, 44.64, 35.89, 25.72))
  expect_equivalent(b0, pdbs$b[, 1:6])

  chain0 <- rbind(rep('A', 6), rep('A', 6),
                  c(NA, rep('A', 5)), c(NA, rep('A', 5)))
  expect_equivalent(chain0, pdbs$chain[, 1:6])

  ali0 <- rbind(c('A', 'R', 'T', 'V', 'K', 'L'), 
                c('A', 'R', 'T', 'V', 'K', 'L'),
                c('-', 'R', 'E', 'V', 'K', 'L'), 
                c('-', 'R', 'E', 'V', 'K', 'L'))
  expect_equivalent(ali0, pdbs$ali[, 1:6])

  sse0 <- rbind(c(' ', ' ', ' ', 'E', 'E', 'E'),
                c(' ', ' ', ' ', 'E', 'E', 'E'),
                c(NA,  ' ', 'E', 'E', 'E', 'E'),
                c(NA,  ' ', 'E', 'E', 'E', 'E'))
  expect_equivalent(sse0, pdbs$sse[, 1:6])

  all.elety0 <- rbind(c('N', 'CA', 'C', 'O', 'CB', 'N'),
                      c('N', 'CA', 'C', 'O', 'CB', 'N'),
                      c(NA,  NA,   NA,  NA,  NA,  'N'),
                      c(NA,  NA,   NA,  NA,  NA,  'N'))
  expect_equivalent(all.elety0, pdbs$all.elety[, 1:6])
 
  all.grpby0 <- c(rep(1, 5), 2)
  expect_equivalent(all.grpby0, pdbs$all.grpby[1:6])

  all.hetatm.xyz0 <- c(8.943,  83.379,  38.955, 2.309,  67.529,  27.056)
  expect_equivalent(all.hetatm.xyz0, pdbs$all.hetatm[[1]]$xyz[, 1:6])

  ## read with multicore
  invisible( capture.output(pdbs.mc <- read.all(aln, ncore=NULL)) )
  expect_equal(pdbs[!names(pdbs) %in% 'call'], pdbs.mc[!names(pdbs.mc) %in% 'call'])

  ## read protein only
  invisible( capture.output(pdbs.prot <- read.all(aln, rm.ligand=TRUE, ncore=NULL)) )
  expect_true(is.null(pdbs.prot$all.hetatm))
  expect_equal(pdbs[!names(pdbs) %in% c('call', 'all.hetatm')],
               pdbs.prot[!names(pdbs.prot) %in% c('call', 'all.hetatm')])
  
  ## read with select atoms
  invisible( capture.output(pdbs.sel <- read.all(aln, sel=c('N', 'CA', 'C'), 
                                             rm.ligand=TRUE, ncore=NULL)) )
  all.elety0 <- rbind(c('N', 'CA', 'C', 'N', 'CA', 'C'),
                      c('N', 'CA', 'C', 'N', 'CA', 'C'),
                      c(NA,  NA,   NA,  'N', 'CA', 'C'),
                      c(NA,  NA,   NA,  'N', 'CA', 'C'))
  expect_equivalent(all.elety0, pdbs.sel$all.elety[, 1:6])

  ## Error handling
  expect_error(read.all('pdb_file_name_is_not_supported'))
  
  ## If one structure is not readable, print warning (ncore=1) and write NAs 
  ## for the missing structure.
  taln <- aln; taln$id[1] <- 'dummy'
  expect_warning(capture.output(pdbs.missing <- read.all(taln, ncore=1)))
  expect_true(all(is.na(pdbs.missing$all[1,])))

  ## If mismatch between structure and alignment, print proper error message
  taln <- aln; taln$ali[1, 16] <- 'A'
  expect_error(capture.output(read.all(taln, ncore=1)))
})


