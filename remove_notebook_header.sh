#!/usr/bin/env bash

# Jupyter is unable to apply the syntax highlighting colors from
# custom themes in jupyterthemes.
# By saving the notebook in the browser with save as the colors
# are preserved but the header remains and the scrolling will be
# disabled.
# These issues can be resolved the way below.

sed -iE 's/<!--\|-->//g;'\
's/\(<div id="header"\)/<!-- \1/;'\
's/\(<div id="site" style=".*">\)/\1 -->/;'\
's/\(id="notebook_panel"\) style="overflow:auto!important;"/\1/;'\
's/\(id="notebook_panel"\)/\1 style="overflow:auto!important;"/'\
 lipyd_demo.html
