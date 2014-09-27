#!/bin/bash

##### SCRIPT TO GENERATE HTML DOCUMENTATIONS WITH STATICDOCS ##
##                                                           ##
## USAGE: ./run_generate_html.sh [yes][no]                   ##
## OPTIONS:                                                  ##
##           yes: Default, execute example codes             ##
##           no:  Don't execute example codes                ##
###############################################################

example=TRUE
if test $# -ge 1; then
   if echo $1 | grep -iq -e false -e no;  then
      example=FALSE
   fi
fi
utildir=`pwd`

# 0. Create a folder for working
workdir=sandbox/`date +%m%d%y.%H%M%S`
mkdir -p $workdir; cd $workdir

# 1. Get clean Bio3D folders
R CMD build $utildir/../bio3d
tar xvfz bio3d*.tar.gz

# 2. Build symbol links; it is necessary when calling devtools:::load_all
ln -s inst/CITATION ./bio3d/
ln -s inst/examples ./bio3d/
ln -s inst/matrices ./bio3d/

# 3. remove dontrun tags to run all example codes
sh $utildir/remove_dontrun.sh

# 4. start an R session and run the commands to generate html files in ./html/
Rscript -e "options(device=x11)" \
        -e "library(staticdocs)" \
        -e "build_site(pkg='bio3d', site_path='html', examples=$example, launch=TRUE)"

# 5. tidy up html files
utildir=$utildir sh $utildir/tidy_html.sh

# 6. refresh your browser

