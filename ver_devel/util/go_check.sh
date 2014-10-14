#!/bin/bash

# Usage: 
#        ./go_check.sh [-cran] [-win] [-nocheck]
#

check=yes
cran=no
win=no
Rscript=Rscript
while test $# -ge 1; do
  if test $1 == "-cran" || test $1 == "-CRAN"; then
     cran=yes
  fi
  if test $1 == "-win" || test $1 == "-WIN"; then
     win=yes
  fi
  if test $1 == "-nocheck"; then
     check=no
  fi
  if test $1 == "-devel"; then
     Rscript=Rscript-devel
  fi
  if test $1 == "-h" || test $1 == "--help"; then
     echo "Usage:"
     echo "       ./go_check.sh [-cran] [-win] [-nocheck]"
     echo 
     exit 0
  fi
  shift
done

if test $check == "yes"; then
   # Use devtools to check as/not as CRAN
   if test $cran == "yes"; then
      # replace devtools::r_env_vars to set NOT_CRAN = "false"
      # it is important to make skip_on_cran() work
      cat > __go_check_tmp.r <<EOF
          r_env_vars <- function() {
              c(R_LIBS = paste(.libPaths(), collapse = .Platform\$path.sep), 
              CYGWIN = "nodosfilewarning", R_TESTS = "", NOT_CRAN = "false", 
              TAR = devtools:::auto_tar())
          }
EOF
      $Rscript -e 'source("__go_check_tmp.r")'\
               -e 'assignInNamespace("r_env_vars", r_env_vars, "devtools")'\
               -e 'devtools::check("../bio3d", cran=TRUE)'
      rm -f __go_check_tmp.r
   else
      $Rscript -e 'devtools::check("../bio3d", cran=FALSE)'
   fi
fi

# Build Windows version and check
# NOTE: it will use winbuilder online and returns results by E-Mail.
if test $win == "yes"; then
   $Rscript -e 'devtools::build_win("../bio3d")'
fi

## Old script to check with default R command s
#R CMD build ../bio3d
#file=(`ls -t bio3d_*.tar.gz`)
#R CMD check --as-cran $file

