#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `lipyd` python module
#
#  Copyright (c) 2015-2018 - EMBL
#
#  File author(s): Dénes Türei (turei.denes@gmail.com)
#
#  This code is not for public use.
#  Please do not redistribute.
#  For permission please contact me.
#
#  Website: http://www.ebi.ac.uk/~denes
#

# from pyteomics import mzml

import pyopenms

class Reader(object):
    
    def __init__(self, fname):
        
        self.fname = fname
        
    def read_mzml(self):
        
        self.fil = pyopenms.MzMLFile()
        self.exp = pyopenms.MSExperiment()
        self.fil.load(self.fname, self.exp)
