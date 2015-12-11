#! /bin/bash

FILES="flens/*.h flens/*.cc flens/*.tcc examples/*.h examples/*.cc"

for i in $FILES; do
    cat $i | ./clean.pl > .tmp
    mv .tmp $i
done
