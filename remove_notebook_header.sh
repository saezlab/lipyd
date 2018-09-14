#!/usr/bin/env bash

sed -iE 's/<!--\|-->//g;'\
's/\(<div id="header"\)/<!-- \1/;'\
's/\(<div id="site" style=".*">\)/\1 -->/;'\
's/\(id="notebook_panel"\) style="overflow:auto!important;"/\1/;'\
's/\(id="notebook_panel"\)/\1 style="overflow:auto!important;"/'\
 lipyd_demo.html
