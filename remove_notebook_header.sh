#!/usr/bin/env bash

sed -iE 's/<!--\|-->//g;s/\(<div id="header"\)/<!-- \1/;s/\(<div id="site" style=".*">\)/\1 -->/' lipyd_demo.html
