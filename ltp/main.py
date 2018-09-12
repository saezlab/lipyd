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
import time

from ltp import table


class ResultsReprocessor(object):
    """
    Reprocesses a directory with top features xlsx tables.
    Finds mgf files for each protein for all protein containing fractions.
    Opens the xlsx files, does custom operations, writes data in new columns
    and saves xlsx files in a new directory with preserving formatting.
    """
    
    def __init__(self, source_dir, mgfdir, target_dir = None):
        
        self.source_dir = source_dir
        self.mgfdir = mgfdir
        self.target_dir = target_dir or '%s__%s' % (
            target_dir,
            time.strftime('%Y.%m.%d_%H.%M')
        )
        
        if not os.path.exists(self.target_dir):
            
            os.mkdir(self.target_dir)
    
    def read_protein_containing_fractions(self):
        
        pass
