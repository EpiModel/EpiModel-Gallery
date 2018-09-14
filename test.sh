#!/bin/bash

d="$(ls -p | grep "/")"

for i in $d; do
  cd $i
  Rscript model.R
  rm *.pdf
  cd ..
done
