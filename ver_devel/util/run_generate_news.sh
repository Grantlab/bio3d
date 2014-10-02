#!/bin/bash

################ SCRIPT TO GENERATE NEWS/CHANGELOG ############
##                                                           ##
## USAGE: ./run_generate_news.sh [start] [end]               ##
## OPTIONS:                                                  ##
##         start: start revision, e.g. v2.0. Default: origin ##
##           end: end revision, e.g. v2.1-0. Default: HEAD   ##
## NOTES:                                                    ##
##      To make this script work, we have to follow          ## 
##      the conventions we made before for commit message,   ## 
##      i.e. four type of commits, "NEW", "ENHANCEMENT",     ##
##      "BUGFIX", and "OTHER", and one empty line between    ##
##      message subject and body.                            ##  
##      See https://bitbucket.org/Grantlab/bio3d/wiki        ##
##      for more details.                                    ##
##                                                           ##
###############################################################

start=
end=HEAD
key=NEW:

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

if test -z $start; then
   range=""
   echo Available versions:
   git tag
else
   range="$start..$end"
fi

git log --no-merges --grep="$key" --diff-filter=A --name-only --date-order \
   --pretty=format:"## %ai | %an | %s" $range ../bio3d/R > log

awk "BEGIN{new=0; print \"$end\\n======\\n\"}
     /^##/{new=1; split(\$0, a, \"|\"); msg=a[length(a)]; sub(/\s*NEW:\s*/, \"\", msg)}
     NF==0{new=0}
     new && !/^##/{split(\$0, f, \"/\"); printf \"* %s%c%s\\n\", f[length(f)], \":\", msg}" \
 log > NEWS

