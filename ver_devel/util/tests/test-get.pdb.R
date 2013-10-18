context("get.pdb(), get pdb file from online")

test_that("get.pdb works properly", {
   ids <- c("1tag", "1tnd", "1tad")  # Gt
   tmp <- tempdir()
   system(paste("rm -f ", tmp, "/*", sep=""))
   expect_identical(get.pdb(ids, tmp, verbose=FALSE), paste(tmp, "/", ids, ".pdb", sep=""))
   expect_warning(get.pdb("3c7kxxx", tmp, verbose=FALSE))
   expect_warning(get.pdb("1tag", tmp, verbose=FALSE))
   expect_identical(get.pdb("1as0", tmp, verbose=FALSE, gzip=TRUE), paste(tmp, "/1as0.pdb", sep=""))
#   expect_error(get.pdb("aaaa", tmp, verbose=FALSE))
   unlink(tmp)
})

test_that("get.pdb, ncore>1 works properly", {
   ids <- c("1tag", "1tnd", "3sn6", "1cul",
            "3v00", "1got", "1kjy", "4g5o", 
            "2ode", "2ihb", "1as0", "1tad")
   tmp1 <- paste(tempdir(), "1", sep="")
   tmp2 <- paste(tempdir(), "2", sep="")
   system(paste("rm -f ", tmp1, "/*", sep=""))
   system(paste("rm -f ", tmp2, "/*", sep=""))
   time1 <- system.time(r1 <- get.pdb(ids, tmp1, verbose=FALSE))
   time2 <- system.time(r2 <- get.pdb(ids, tmp2, ncore=2, verbose=FALSE))
   expect_identical(r2, paste(tmp2, "/", ids, ".pdb", sep=""))
   expect_identical(list.files(tmp1), list.files(tmp2))
   if("multicore" %in% installed.packages() & 
     multicore:::detectCores() > 1)
      expect_true(time1["elapsed"] >= time2["elapsed"])
   unlink(tmp1)
   unlink(tmp2)
})

