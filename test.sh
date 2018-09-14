#!/bin/bash

set -e

d="$(ls -p | grep "/")"

for i in $d; do
  cd $i
  Rscript model.R options(error = function() q("no", 1, FALSE))
  rm *.pdf
  cd ..
done
