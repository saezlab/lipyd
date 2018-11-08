#!/usr/bin/bash

BASEDIR='/home/denes/archive/ltp'

find $BASEDIR -depth -not -path '*deleteme*' -path '*/Results/*mgf' -exec ln -s -f {} $BASEDIR/mgf_invivo/ \;
find $BASEDIR -depth -not -path '*deleteme*' -path '*update*/Results/*mgf' -exec ln -s -f {} $BASEDIR/mgf_invivo/ \;

# fix file name typos
for f in `ls $BASEDIR/*ORP9STARD15*`; do mv $f `echo $f | sed 's/STARD15//g'`; done
