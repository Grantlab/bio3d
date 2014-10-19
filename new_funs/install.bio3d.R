
install.bio3d <- function(version="2.1-0", dependencies=TRUE, cran.repos="http://cran.at.r-project.org/") {
  
  ## ToDo:
  ##      - Check if we have permissions to install.
  ##      - Use 'try()' to catch warning messages for loading packages that are not yet installed.
  ##      - Test robustness of checking for muscle and dssp on windows.
  ##      - Add windows specific binary package installation (currently its source only) 
  ##      - Implement an update option that checks on version installed vs that requested.
  ##         Currently this function will install the version requested over whatever is already there with no checking.

  ##- Check if we have permissions to install
  # To Do!

  ##- Check and install missing R package dependences if necessary
  if(dependencies){
    packages <- c("ncdf", "lattice", "bigmemory", "XML", "RCurl", "igraph", "knitr")

    get.package <- function(x, cran.repos=cran.repos){
      if (!require(x, character.only = TRUE)){
        install.packages(x,, dependencies=TRUE, repos=cran.repos)
          if(!require(x, character.only = TRUE)){
            stop(paste("Dependent package",x,"not found on CRAN:",cran.repos))
          }
      } else {
      cat("Dependent package",x,"already installed\n")
      }
    }
    packages.rtn <- sapply(packages, get.package, cran.repos)
  }


  ##- Download and install Bio3D from the package website
  online.file <- paste0("http://thegrantlab.org/bio3d/phocadownload/Bio3D_version",
    substr(version,1,1),".x/bio3d_",version,".tar.gz")
  local.file <- basename(online.file)
  status <- download.file(online.file, destfile=local.file)
  install.packages(local.file, repos=NULL, type="source")
  unlink(local.file)


  ##- Check on missing utility programs (dssp and muscle)
  utilities <- c("dssp", "muscle")#, "stride", "mustang", "makeup")
  missing.util <- nchar(Sys.which(utilities)) == 0
  if( any(missing.util) ) {
    warning(paste0("  Checking for external utility programs failed\n",
      "    Please make sure '", paste(names(missing.util[missing.util]), collapse="', '"),
      "' is in your search path, see:\n",
      "    http://thegrantlab.org/bio3d/tutorials/installing-bio3d#utilities"))
  } else {
    cat("External utility programs found\n")
  }

  ##- Alternate check approach - not yet sure of the best way to do this robustly
  #  os1 <- .Platform$OS.type
  #  if (os1 == "windows") {
  #    status.muscle <- shell(shQuote("muscle -version"), ignore.stderr=TRUE, ignore.stdout=TRUE)
  #    status.dssp <- shell(shQuote("dssp --version"), ignore.stderr=TRUE, ignore.stdout=TRUE)
  #  } else {
  #    status.muscle <- system("muscle -version", ignore.stderr = TRUE, ignore.stdout = TRUE)
  #    status.dssp <- system("dssp --version",ignore.stderr = TRUE, ignore.stdout = TRUE)
  #  }
  #  if (!(status.muscle %in% c(0, 1)))
  #      warning(paste0("Checking for external program failed\n", 
  #        "  make sure 'muscle' is in your search path, see: http://thegrantlab.org/bio3d/tutorials/installing-bio3d#utilities") )
  #
  #  if (!(status.dssp %in% c(0, 1)))
  #      warning(paste0("Checking for external program failed\n", 
  #        "  make sure 'dssp' is in your search path, see http://thegrantlab.org/bio3d/tutorials/installing-bio3d#utilities"))
  #}


}


