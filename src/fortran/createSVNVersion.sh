#!/bin/bash

# create a sed script containing svn version

OUT=svnversion.sed

echo "s|svn_vers_string|$(svnversion)|" > $OUT.tmp

if [ -f $OUT ]; then
  if [ -n "$(diff --brief $OUT $OUT.tmp)" ]; then
     cp $OUT.tmp $OUT
  fi
else
  cp $OUT.tmp $OUT
fi

rm -f $OUT.tmp
