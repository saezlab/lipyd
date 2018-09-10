#!/usr/bin/bash

find . -depth -not -path '*deleteme*' -path '*/Results/*mgf' -exec ln -s {} ./mgf_invivo/ \;
find . -depth -not -path '*deleteme*' -path '*update*/Results/*mgf' -exec ln -s -f {} ./mgf_invivo/ \;

# fix file name typos
for f in `ls *ORP9STARD15*`; do mv $f `echo $f | sed 's/STARD15//g'`; done
