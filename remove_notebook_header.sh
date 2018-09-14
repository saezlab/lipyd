!#/usr/bin/env bash

sed -ire 's/\(<div id="header"\)/<!-- \1/;s/\(<div id="site" style=".*">\)/\1 -->/' lipyd_demo.html
