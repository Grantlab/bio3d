#!/bin/bash

################ SCRIPT TO GENERATE NEWS/CHANGELOG ############
##                                                           ##
## USAGE: ./run_generate_news.sh [start] [end] [key]         ##
## OPTIONS:                                                  ##
##         start: start revision, e.g. v2.0. Default: origin ##
##           end: end revision, e.g. v2.1-0. Default: HEAD   ##
##           key: one of the four types of commits.          ##
##                By default, NEW.                           ##
## NOTES:                                                    ##
##      To make this script work, we have to follow          ## 
##      the conventions we made before for commit message,   ## 
##      i.e. four type of commits, "NEW", "ENHANCEMENT",     ##
##      "BUGFIX", and "OTHER", and two empty lines between   ##
##      message subject and body.                            ##  
##      See https://bitbucket.org/Grantlab/bio3d/wiki        ##
##      for more details.                                    ##
##                                                           ##
###############################################################

start=
end=
key=NEW

if test $# -gt 0; then 
   start=$1
   shift
fi
if test $# -gt 0; then 
   end=$1
   shift
fi

if test $# -gt 0; then 
   key=$1
   shift
fi

git log --no-merges --grep="$key:" --stat --pretty=format:"%an: %s%n%b" $start..$end > NEWS

