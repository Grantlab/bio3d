#!/bin/bash


## This script generates a new directory in which should be copied to the 
## Joomla! website. It substitutes the path of the figures in the .md file 
## to a path convienent for use on the web-site. It also copies
## the PNG files in the figure-markdown directory to 'figures'. 


if test $# -lt 1; then
   echo Usage: ./markdown_for_web.sh prefix-of-md-file
   exit 0
fi


prefix=$1
file=${prefix}".md"
fromdir=${prefix}"_files\/figure-markdown_phpextra"
fromdir2=${prefix}"_files/figure-markdown_phpextra"
finaldir=${prefix}
figdir="figures"

echo " "

if [ ! -s $file ]; then
    echo "File not found: "${file}
    echo Usage: ./markdown_for_web.sh prefix-of-md-file
    exit 0;
fi

if [ ! -d ${fromdir2} ]; then
    echo "Directory not found: "${fromdir2}
    exit 0;
fi

if [ -d ${finaldir} ]; then
    echo "... Output directory '$finaldir' already exists"
    exit 0;
fi

if [ ! -d ${figdir} ]; then
    echo "... NOTE: Figures directory ${figdir} does not exist"
fi

if [ ! -d ${finaldir} ]; then
    echo "... Making directory '$finaldir'"
    mkdir ${finaldir}
fi

str="s/${fromdir}/${figdir}/g"


echo "... Running command: \"sed '$str' $file > ${finaldir}/${file} \""
sed $str $file > ${finaldir}"/"${file}

if [ ! -s ${finaldir}"/"${file} ]; then
    echo "... Error: A problem occured with sed command"
    exit 0;
fi

if [ ! -d ${finaldir}/${figdir} ]; then
    echo "... Making directory '${finaldir}/${figdir}'"
    mkdir ${finaldir}/${figdir}
fi

if [ -d ${figdir} ]; then
    echo "... Copying figure PNG files from '$figdir' to '${finaldir}/${figdir}'"
    cp -p ${figdir}/*.png ${finaldir}/${figdir}
fi

if [ -d ${fromdir2} ]; then
    echo "... Copying figure PNG files from '$fromdir2' to '${finaldir}/${figdir}/'"
    cp -p ${fromdir2}/*.png ${finaldir}/${figdir}/
fi


echo " "
echo "========================== "
echo " "
echo "Copy the following directory to the web-server:"
echo "- "${finaldir}"/"
echo " "
echo "========================== "
echo " "

