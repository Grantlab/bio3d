check.bio3d <- function(pkg, as.cran = TRUE, run.dontrun = FALSE, run.donttest =FALSE,
   cran.repos = "http://cran.at.r-project.org/", cleanup = FALSE, ...) {

   ##- Check if pkg is a file or a directory
   if(!is.character(pkg))
      stop("Provide path to bio3d source file or directory")
   if(length(list.files(pkg)) == 0) {
      tp = tempdir()
      unlink(file.path(tp, "bio3d"), recursive = TRUE) 
      if(untar(pkg, exdir = tp)) {
         stop("Unable to extract package source from file")
      }
      pkg = file.path(tp, "bio3d")
   }
  
   ##- Check and install dependencies
   packages <- c("devtools", "roxygen2", "ncdf", "lattice", 
                 "XML", "RCurl", "igraph", "knitr", "testthat")
   if(.Platform$OS.type == "unix") packages <- c(packages, "bigmemory")
   get.package <- function(x, cran.repos = cran.repos){
     if (!requireNamespace(x, quietly = TRUE)){
        install.packages(x, dependencies = TRUE, repos = cran.repos)
        if(!requireNamespace(x, quietly = TRUE)){
           stop(paste("Dependent package", x, "not found on CRAN:", cran.repos))
        }
     } else {
       cat("Dependent package", x, "already installed\n")
     }
   }
   packages.rtn <- sapply(packages, get.package, cran.repos)
 
   if(as.cran) {
      # replace devtools::r_env_vars to set NOT_CRAN = "false"
      # it is important to make skip_on_cran() work
      r_env_vars <- function() {
         c(R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep), 
         CYGWIN = "nodosfilewarning", R_TESTS = "", NOT_CRAN = "false", 
         TAR = devtools:::auto_tar())
      }
      assignInNamespace("r_env_vars", r_env_vars, "devtools")
   }
   
   args = NULL 
   if(run.dontrun) args <- c(args, "--run-dontrun")
   if(run.donttest) args <- c(args, "--run-donttest")

   devtools::check(pkg, cran=as.cran, document=FALSE, force_suggests=FALSE, 
      args = args, cleanup = cleanup, ...) 
}
