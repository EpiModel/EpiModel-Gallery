#!/bin/bash

set -e

d="$(ls -p | grep "/")"
dl="$(find ./* -maxdepth 0 -type d | wc -l)"

echo " "
echo "#############################################"
echo "EpiModel Gallery Testing:" $dl "Directories"
echo "---------------------------------------------"

for i in $d; do
  cd $i
  echo -n $i "... "
  SECONDS=0
  Rscript model.R "options(error = function() q('no', 1, FALSE))" >& /dev/null
  echo -n "OK" "... Elapsed Time: $SECONDS seconds."
  rm *.pdf
  cd ..
done

echo "#############################################"
echo " "
