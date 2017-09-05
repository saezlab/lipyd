#!/bin/bash

rm src/*.pyc
rm src/emese/*.pyc
rm -r src/emese.egg-info

sed -i 's/\([0-9]*\.[0-9]*\.\)\([0-9]*\)/echo \1$\(\(\2+1\)\)/ge' src/emese/__version__

python setup.py sdist
