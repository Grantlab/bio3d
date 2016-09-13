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

# 1. Create a folder for working
workdir=sandbox/`date +%m%d%y.%H%M%S`
mkdir -p $workdir; cd $workdir

# 2. Get clean Bio3D folders
R CMD build $utildir/../bio3d
tar xvfz bio3d*.tar.gz

# 3. Build symbol links; it is necessary when calling devtools:::load_all
ln -s inst/CITATION ./bio3d/
ln -s inst/examples ./bio3d/
ln -s inst/matrices ./bio3d/

# 4. remove dontrun tags to run all example codes
sh $utildir/remove_dontrun.sh

# 5. start an R session and run the commands to generate html files in ./html/
mkdir html
if ! Rscript -e "library(staticdocs)" \
        -e "options(device=x11)" \
        -e "build_site(pkg='bio3d', site_path='html', examples=$example, launch=TRUE)"; then
   echo "Error: running staticdocs"
   exit 1
fi

# 6. tidy up html files
utildir=$utildir sh $utildir/tidy_html.sh
# some modification to index.html...
# - add links to all vignettes
sed -i 's/<li><a href="vignettes\/bio3d_vignettes.html">bio3d Vignettes<\/a><\/li>/<li><a href="vignettes\/Bio3D_install.pdf">Installing Bio3D<\/a><\/li>\
      <li><a href="vignettes\/Bio3D_pdb.pdf">PDB Structure Manipulation and Analysis<\/a><\/li>\
      <li><a href="vignettes\/Bio3D_pca.pdf">Comparative Sequence and Structure Analysis with Bio3D<\/a><\/li>\
      <li><a href="vignettes\/Bio3D_nma.pdf">Enhanced Methods for Normal Mode Analysis with Bio3D<\/a><\/li>\
      <li><a href="vignettes\/Bio3D_md.pdf">Beginning Trajectory Analysis with Bio3D<\/a><\/li>\
      <li><a href="vignettes\/cna_vignette.pdf">Protein Structure Networks with Bio3D<\/a><\/li>/' html/index.html 

# 7. create a link to the results
rm -f $utildir/html
ln -s $workdir/html $utildir

# 8. refresh your browser

