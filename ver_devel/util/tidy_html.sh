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

## if two figures placed on continuous lines, delete the first one
#for i in html/*.html; do 
##   echo $i
#   awk 'BEGINE{n=-1} NR>1 && (!/png/||NR>n+1) {print a} \
#     /png/{n=NR} {a=$0} END{print a}' $i \
#    > `basename $i`.t
#   mv `basename $i`.t $i
#done

#./check_html.sh > tidy_html.log

# Remove additional vignettes; -- expired with new version of staticdocs!
#lines=(`awk '/href=\"vignettes\//{printf "%d ", NR}' html/index.html`)
#nums=`expr ${#lines[@]} / 2`
#for((i=0; i<$nums; i++)); do
#   lines=(`awk '/href=\"vignettes\//{printf "%d ", NR}' html/index.html`)
#   sed "${lines[0]}"d html/index.html > index.html
#   mv index.html html/index.html
#done

#cp $utildir/../bio3d/vignettes/*.pdf html/vignettes/
# copy link to website
rm -rf html/vignettes
ln -s ../phocadownload/vignettes html/vignettes

# Remove multiple lines of progress bar
for i in html/*.html; do
   sed 's/.*\(\s*|=*|\s*100%\)/\n\1/' $i > t.html
   mv t.html $i
done

# remove returns in the index.html
sed 's/<br \/>/:/' html/index.html > index.html
mv index.html html/

# additional output file
if test -r eg.html; then  mv eg.html html/; fi

# remove messy "skipping download" warning
for fnm in `grep -in 'Skipping download' html/*.html | awk '{split($0, a, ":"); print a[1]}'`; do
   l2=$( grep -in -m 1 'Skipping download' $fnm | cut -f1 -d: )
   let l1=l2-1
   echo Deleting lines: $fnm, $l1, $l2
   sed -i "$l1,$l2"d $fnm
done


# remove very messy "cannot wait for child ..." warning because of older version of parallel
for fnm in `grep -in 'cannot wait for child' html/*.html | awk '{split($0, a, ":"); print a[1]}'`; do
   l2=$( grep -in -m 1 'cannot wait for child' $fnm | cut -f1 -d: )
   let l1=l2-1
   echo Deleting lines: $fnm, $l1, $l2
   sed -i "$l1,$l2"d $fnm
done
