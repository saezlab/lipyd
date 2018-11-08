#!/bin/bash

rm src/*.pyc
rm src/lipyd/*.pyc
rm -r src/lipyd.egg-info

sed -i 's/\([0-9]*\.[0-9]*\.\)\([0-9]*\)/echo \1$\(\(\2+1\)\)/ge' src/lipyd/__version__

python setup.py sdist
