#!/bin/bash
##### This script removes the dontrun tags and so
##### all the example codes will be executed

for i in ../bio3d/man/*.Rd; do
   # skip some examples because of the missing data, run errors
   # or just to avoid confusion of using attach(data) (pca.xyz.Rd)
   # we may think of another consistent example in future 
   if [ `basename $i` != "rmsip.Rd" ] && \
      [ `basename $i` != "pca.xyz.Rd" ] && \
      [ `basename $i` != "read.mol2.Rd" ] && \
      [ `basename $i` != "bio3d.package.Rd" ] && \
      [ `basename $i` != "binding.site.Rd" ] && \
      [ `basename $i` != "plot.bio3d.Rd" ]; then
      echo $i
      # find and delete tags \dontrun{ and }
      awk 'BEGIN{bok=0; n=0} 
           /\\dontrun\{/{bok=1}
           bok {a=$0; while(sub(/\{/,"",a)) n++; while(sub(/\}/,"",a)) n--}
           n>0 && !/\\dontrun\{/ {print} 
           !bok{print} 
           bok && n==0 {bok=0}' $i > t.Rd
      mv t.Rd $i
   fi
done
