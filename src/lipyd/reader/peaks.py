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

import csv

import reader.xls


class PeaksReader(object):
    
    def __init__(self, fname, ionmode = None, format = None):
        
        self.fname = fname
        
        self.guess_ionmode(ionmode)
        
    
    def read(self):
        
        with open(self.fname, 'r') as fp:
            
            hdr = fp.readline()
