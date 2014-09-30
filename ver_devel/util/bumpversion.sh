#!/bin/bash

if test $# -lt 1; then
   echo Usage: ./bumpversion.sh new-version
   exit 0
fi

ver=$1

# check branch
if test `git br | awk '/^*/{print $2}'` != "release"; then
   echo Warning: You are not on \"release\" branch
   echo Switch automatically...
   if ! git co release; then
      echo Failed to switch...stop
      exit 1
   fi
fi

# update versions
sed -i "s/\(^Version:\s*\).*/\1$ver/" ../bio3d/DESCRIPTION 
sed -i "s/\(^Version:\s*\\\tab\s*\).*\(\s*\\\cr\)/\1$ver\2/" ../bio3d/man/bio3d.package.Rd
sed -i "s/\(^Date:\s*\\\tab\s*\).*\(\s*\\\cr\)/\1`date +'%Y-%m-%d'`\2/" ../bio3d/man/bio3d.package.Rd

