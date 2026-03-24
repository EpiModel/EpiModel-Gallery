#!/bin/bash

set -e

d="$(ls -p examples/ | grep "/")"
dl="$(echo "$d" | wc -l | tr -d ' ')"

echo " "
echo "#############################################"
echo "EpiModel Gallery Testing:" $dl "Examples"
echo "---------------------------------------------"

for i in $d; do
    dir="examples/$i"
    echo -n "$dir ... "
    SECONDS=0
    if Rscript "${dir}model.R" "options(error = function() q('no', 1, FALSE))" >& /dev/null
    then
        echo "OK" "($SECONDS seconds)"
    else
        test_failed=1
        echo "Failed"
    fi
    cd "$dir"
    rm -f *.pdf
    cd ../..
done

rm -f *.pdf

echo "#############################################"
echo " "

if [ -z ${test_failed+x} ]
then
    echo "All Tests Successful"
    exit 0
else
    echo "Some Tests Failed"
    exit 1
fi
