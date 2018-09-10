#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `ltp` python module
#
#  Copyright (c) 2018 - EMBL
#
#  File author(s): Dénes Türei (turei.denes@gmail.com)
#
#  This code is not for public use.
#  Please do not redistribute.
#  For permission please contact me.
#
#  Website: http://www.ebi.ac.uk/~denes
#

import os

from ltp import table

class ReadResults(object):
    
    def __init__(self, result_dir, mgfdir):
        
        self.result_dir = result_dir
        self.mgfdir = mgfdir
