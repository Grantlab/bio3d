#!/bin/bash

# remove 0byte figures
empfiles=(`ls -l html/*.png | awk '$5==0{print $9}'`)
if [ ${#empfiles[@]} -gt 0 ]; then
   for i in ${empfiles[@]}; do
      for j in html/*.html; do 
         sed "/`basename $i`/"d $j > `basename $j`.t 
         mv `basename $j`.t $j
#         rm $i
      done 
   done
fi

# remove figures doesn't exist 
for j in html/*.html; do 
   pngfiles=(`grep '<img src=.*png' $j | sed "s/^.*src='//" | sed "s/'.*//"`)
   if [ ${#pngfiles[@]} -gt 0 ]; then
      for i in ${pngfiles[@]}; do
         if [ ! -r html/$i ]; then
            sed "/`basename $i`/"d $j > `basename $j`.t 
            mv `basename $j`.t $j
         fi
      done 
   fi
done

# if two figures placed on continuous lines, delete the first one
for i in html/*.html; do 
#   echo $i
   awk 'BEGINE{n=-1} NR>1 && (!/png/||NR>n+1) {print a} \
     /png/{n=NR} {a=$0} END{print a}' $i \
    > `basename $i`.t
   mv `basename $i`.t $i
done

#./check_html.sh > tidy_html.log

# remove return in the index.html
sed 's/<br \/>/:/' html/index.html > index.html
mv index.html html/

mv eg.html html/

