#!/bin/bash
##### This script removes the dontrun tags and so
##### all the example codes will be executed

for i in ./bio3d/man/*.Rd; do
   # skip some examples because of the missing 
   # data or run errors
#   if [ `basename $i` != "read.crd.Rd" ] && \
#      [ `basename $i` != "read.crd.amber.Rd" ] && \
#      [ `basename $i` != "read.crd.charmm.Rd" ] && \
#      [ `basename $i` != "as.pdb.Rd" ] && \
#      [ `basename $i` != "atom.select.Rd" ] && \
#      [ `basename $i` != "read.prmtop.Rd" ] && \
#      [ `basename $i` != "read.mol2.Rd" ]; then
      echo $i
      # find and delete tags \dontrun{ and }
      awk 'BEGIN{bok=0; n=0} 
           /\\dontrun\{/{bok=1}
           bok {a=$0; while(sub(/\{/,"",a)) n++; while(sub(/\}/,"",a)) n--}
           n>0 && !/\\dontrun\{/ {print} 
           !bok{print} 
           bok && n==0 {bok=0}' $i > t.Rd
      mv t.Rd $i
#   fi
    
   # find help(), demo(), identify or identify.cna() and add \dontrun{}
   sed -e '/\\examples\s*{/,$s/^\(\s*help\s*(.*).*\)$/\\dontrun{\n\1\n}/' \
       -e '/\\examples\s*{/,$s/^\(\s*demo\s*(.*).*\)$/\\dontrun{\n\1\n}/' \
       -e '/\\examples\s*{/,$s/^\([^#]*identify\s*(.*).*\)$/\\dontrun{\n\1\n}/' \
       -e '/\\examples\s*{/,$s/^\([^#]*identify.cna\s*(.*).*\)$/\\dontrun{\n\1\n}/' \
   $i > t.Rd
   mv t.Rd $i

   if [ `basename $i` == "hmmer.Rd" ]; then
     # replace the tag '\cr' to newlines to avoid some errors from staticdocs
     sed -i 's/\\cr//g' $i
   fi

   # comment out 'check.utility()' in examples, because it will cause
   #  the loss of some plots.
   if [ `basename $i` != "check.utility.Rd" ]; then
      awk 'BEGIN{bok=0; n=0; bok_inv=0} 
           /check.utility/{bok=1}
           /\!check.utility/{bok_inv=1}
           bok {a=$0; while(sub(/\{/,"",a)) n++; while(sub(/\}/,"",a)) n--}
           !bok_inv && n>0 && !/check.utility/ && !(/else/ && n==1){print} 
           /else/ && n==1 {bok_inv=0}
           !bok{print} 
           bok && n==0 {bok=0}' $i > t.Rd
      mv t.Rd $i
   fi
   
   # comment out 'requireNamespace()' in examples, because it will cause
   #  the loss of some plots.
#   if [ `basename $i` != "check.utility.Rd" ]; then
      awk 'BEGIN{bok=0; n=0; bok_inv=0}
           /requireNamespace/{bok=1}
           /\!requireNamespace/{bok_inv=1}
           bok {a=$0; while(sub(/\{/,"",a)) n++; while(sub(/\}/,"",a)) n--}
           !bok_inv && n>0 && !/requireNamespace/ && !(/else/ && n==1){print} 
           /else/ && n==1 {bok_inv=0}
           !bok{print} 
           bok && n==0 {bok=0}' $i > t.Rd
      mv t.Rd $i
#   fi
done
