#!/bin/bash

set -e

d="$(ls -p | grep "/")"

echo " "
echo "#############################################"
echo "Starting EpiModel Gallery Testing"
echo "---------------------------------------------"

for i in $d; do
  cd $i
  echo -n $i "... "
  Rscript model.R "options(error = function() q('no', 1, FALSE))" >& /dev/null
  echo "OK"
  rm *.pdf
  cd ..
done

echo "#############################################"
echo " "
