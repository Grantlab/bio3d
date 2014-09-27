#!/bin/bash

##### STEPS TO GENERATE HTML DOCS WITH STATICDOCS #######
##### MUST BE RUN UNDER VER_DEVEL/UTIL/ !!!

if [ $# -lt 1 ]; then
   echo THIS SCRIPT NEED TO RUN INTERACTIVELY
   exit 1
fi

utildir=`pwd`

# 0. Create a folder for working
workdir=sandbox/`date +%m%d%y.%H%M%S`
mkdir -p $workdir; cd $workdir

# 1. Backup docs
R CMD build $utildir/../bio3d
tar xvfz bio3d*.tar.gz

# 2. Build symbol links; it is necessary when calling devtools:::load_all
ln -s inst/CITATION ./bio3d/
ln -s inst/examples ./bio3d/
ln -s inst/matrices ./bio3d/

# 3. remove dontrun tags to run all example codes
sh $utildir/run_remove_dontrun.sh

# 4. run staticdocs and produce html files in html/
utildir=$utildir R
#    The script must be run in R interactively
source(file.path(Sys.getenv("utildir"), "run_staticdocs.r"))

# 5. tidy up html files
./run_tidy_html.sh

# 6. Restore docs
#rm -rf ./bio3d/man
#mv ./man_bak ./bio3d/man

# 7. Remove the symbol links
#rm ./bio3d/CITATION
#rm ./bio3d/examples
#rm ./bio3d/matrices
