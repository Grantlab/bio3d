
install.bio3d <- function(version="2.1-0", dependencies=TRUE, cran.repos="http://cran.at.r-project.org/") {
  
  ## ToDo: implement checking for muscle and dssp.
  ##       implement an update option.

  ##- Check and install missing R package dependences if necessary
  if(dependencies){
    packages <- c("ncdf", "lattice", "grid", 
                  "bigmemory", "parallel", "XML",
                  "RCurl", "igraph")

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

  ##- Download and install Bio3D from package website
  online.file <- paste0("http://thegrantlab.org/bio3d/phocadownload/Bio3D_version",
    substr(version,1,1),".x/bio3d_",version,".tar.gz")
  local.file <- basename(online.file)
  status <- download.file(online.file, destfile=local.file)
  install.packages(local.file, repos=NULL, type="source")
  unlink(local.file)

  ##- Check for external working DSSP and muscle
  if(FALSE){ ##-- Not yet sure of the best way to do this robustly ...
    os1 <- .Platform$OS.type
    #if (os1 == "windows") {
    #  status.muscle <- shell(shQuote("muscle -version"), ignore.stderr=TRUE, ignore.stdout=TRUE)
    #  status.dssp <- shell(shQuote("dssp --version"), ignore.stderr=TRUE, ignore.stdout=TRUE)
    #} else {
    status.muscle <- system("muscle -version", ignore.stderr = TRUE, ignore.stdout = TRUE)
    status.dssp <- system("dssp --version",ignore.stderr = TRUE, ignore.stdout = TRUE)
    #}
    if (!(status.muscle %in% c(0, 1)))
        warning(paste0("Checking for external program failed\n", 
          "  make sure 'muscle' is in your search path, see: http://thegrantlab.org/bio3d/tutorials/installing-bio3d") )

    if (!(status.dssp %in% c(0, 1)))
        warning(paste0("Checking for external program failed\n", 
          "  make sure 'dssp' is in your search path, see http://thegrantlab.org/bio3d/tutorials/installing-bio3d"))
  }
}


