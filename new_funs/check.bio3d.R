check.bio3d <- function(pkg, run.skipcran = FALSE, run.dontrun = FALSE, run.donttest =FALSE,
   example.only = FALSE, start = NULL, cran.repos = "http://cran.rstudio.com/", 
   cleanup = FALSE, ...) {

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
 
   ##- Install devtools
   if(!requireNamespace('devtools', quietly = TRUE)) {
      install.packages('devtools', dependencies = TRUE, repos = cran.repos)
   } else {
       cat("devtools already installed\n")
   }

   ##- Check and install dependencies
   devtools::install_deps(pkg, dependencies = TRUE)

   if(example.only) {
      devtools::run_examples(pkg, start = start, show = TRUE, test = !run.donttest, 
          run = !run.dontrun, fresh = FALSE)
   } else {
      if(!run.skipcran) {
         # replace devtools::r_env_vars to set NOT_CRAN = "false"
         # it is important to make skip_on_cran() work
         r_env_vars0 <- devtools::r_env_vars
         r_env_vars <- function() {
            c(R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep), 
            CYGWIN = "nodosfilewarning", R_TESTS = "", NOT_CRAN = "false", 
            TAR = devtools:::auto_tar())
         }
         assignInNamespace("r_env_vars", r_env_vars, "devtools")
         on.exit(assignInNamespace("r_env_vars", r_env_vars0, "devtools")) 
      }
      
      args = NULL 
      if(run.dontrun) args <- c(args, "--run-dontrun")
      if(run.donttest) args <- c(args, "--run-donttest")
   
      devtools::check(pkg, cran=TRUE, document=FALSE, force_suggests=FALSE, 
         args = args, cleanup = cleanup, ...)
   }  
}
