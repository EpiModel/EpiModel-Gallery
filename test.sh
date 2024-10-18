#!/bin/bash

set -e

d="$(ls -p | grep "/")"
dl="$(find ./* -maxdepth 0 -type d | wc -l)"

echo " "
echo "#############################################"
echo "EpiModel Gallery Testing:" $dl "Directories"
echo "---------------------------------------------"

for i in $d; do
  echo -n $i "... "
  SECONDS=0
  if Rscript "$i/model.R" "options(error = function() q('no', 1, FALSE))" >& /dev/null
  then
    echo "OK" "($SECONDS seconds)"
  else
    echo "Failed"
  fi
done

echo "#############################################"
echo " "
