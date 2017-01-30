#!/usr/bin/env python2
# -*- coding: utf-8 -*-

#
#  Copyright (c) 2015-2017 - EMBL-EBI
#
#  File author(s): Dénes Türei (denes@ebi.ac.uk)
#
#  This code is not for public use.
#  Please do not redistribute.
#  For permission please contact me.
#
#  Website: http://www.ebi.ac.uk/~denes
#

#
    # filtering features accross samples and fractions
    # to find those relevant ones 
    # belonging to LTP bound lipids
#

import numpy as np
import scipy as sp
import time
import os
import sys
import copy
import pickle

from future.utils import iteritems

#from common import *

#

import emese

#l = ltp2.LTP(recal_source = 'denes', only_marcos_fragments = False, marco_lipnames_from_db = False)
l = emese.Screening(only_marcos_fragments = False)


'''
Initializing from the beginning.
'''
l.init_from_scratch()

'''
Here we reload the workspace including all the data from the pickle:
'''
#l.init_reinit()

#l.load_data()

'''
Apply basic filters, and obtain the valid features.
'''
l.valid_features(cache = False)
'''
Processing standards from mzML
Calculating drifts
Recalibrating all m/z's
'''
#l.recalibration()

'''
Lookup lipids for all the valid features. (variable name: `lip`)
'''
l.ms1()
l.ms2_onebyone()
