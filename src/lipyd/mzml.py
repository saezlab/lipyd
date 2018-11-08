#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `lipyd` python module
#
#  Copyright (c) 2015-2018 - EMBL
#
#  File author(s): Dénes Türei (turei.denes@gmail.com)
#
#  Distributed under the GNU GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://www.ebi.ac.uk/~denes
#

# from pyteomics import mzml

try:
    import pyopenms
except:
    pass

class Reader(object):
    
    def __init__(self, fname):
        
        self.fname = fname
        
    def read_mzml(self):
        
        self.fil = pyopenms.MzMLFile()
        self.exp = pyopenms.MSExperiment()
        self.fil.load(self.fname, self.exp)
