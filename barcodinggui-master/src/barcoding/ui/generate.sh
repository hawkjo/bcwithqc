#!/bin/sh

UIC="uic-qt6"
if ! which $UIC > /dev/null 2>&1 ; then
    UIC="/usr/lib/qt6/uic"
fi
for f in *.ui ; do
    $UIC -g python -o `basename $f .ui`.py $f
done
