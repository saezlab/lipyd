#!/usr/bin/env python
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
import cProfile as profile
import pyprof2calltree

from future.utils import iteritems

#from common import *

#

import emese

#l = ltp2.LTP(recal_source = 'denes', only_marcos_fragments = False, marco_lipnames_from_db = False)
l = emese.Screening(only_marcos_fragments = True,
                    slope_profile_selection = True,
                    pcont_fracs_from_abs = True)

l = emese.Screening(slope_profile_selection = False,
                    sec_filenames_protein_first = False,
                    filter_wrong_fractions = False,
                    fractionsf = 'control_sample.csv',
                    fr_offsets = [0.010, 0.045],
                    use_gel_profiles = False,
                    pp_do_correction = True,
                    fixed_fraction_layout = True)
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

p2 = profile.Profile()
p2.enable()
l.ms2_onebyone()
p2.disable()
p2.dump_stats('ms2_profiler_20170301')
pyprof2calltree.visualize(p2.getstats())

# profiling MS2
import cProfile as profile

p = profile.Profile()
p.enable()
l.ms2_oneprotein('STARD1')
p.disable()

l.read_manual()
l.manual_df(screen = 'E')
l.export_df(l.pmanual, 'enric_processed.csv')

l.read_manual2()
l.manual_df(screen = 'A', only_swl_col = True)
l.export_df(l.pmanual, 'antonella_final.csv')

l.combined_network()
l.combined_network(classes = set(['I']), fname = 'combined_network_classI_%s.csv')


l.fractions_barplot2(features = True,
                     highlight = 'enrbb',
                     ignore_offsets = True,
                     all_features = False,
                     pdfname = 'pp_features_enric.pdf')

l.peak_ratio_score_threshold = 1.5
l.peak_ratio_score_bool()
l.fractions_barplot2(features = True,
                     highlight = 'prs1',
                     ignore_offsets = True,
                     all_features = False,
                     pdfname = 'pp_features_prs0.5.pdf')

l.logical_and(one = 'prs1', two = 'enrbb', result = 'prs05_enr')
l.fractions_barplot2(features = True,
                     highlight = 'prs05_enr',
                     ignore_offsets = True,
                     all_features = False,
                     pdfname = 'pp_features_enric_prs0.5.pdf')

l.peak_ratio_score_threshold = 1.5
l.peak_ratio_score_bool()
l.logical_and(one = 'prs1', two = 'enrbb', result = 'prs15_enr')
l.fractions_barplot2(features = True,
                     highlight = 'prs15_enr',
                     ignore_offsets = True,
                     all_features = False,
                     pdfname = 'pp_features_enric_prs1.5.pdf')

l.peak_ratio_score_threshold = 1.0
l.peak_ratio_score_bool()
l.logical_and(one = 'prs1', two = 'enrbb', result = 'prs10_enr')
l.fractions_barplot2(features = True,
                     highlight = 'prs10_enr',
                     ignore_offsets = True,
                     all_features = False,
                     pdfname = 'pp_features_enric_prs1.0.pdf')

l.peak_ratio_score_threshold = 0.7
l.peak_ratio_score_bool()
l.logical_and(one = 'prs1', two = 'enrbb', result = 'prs07_enr')
l.fractions_barplot2(features = True,
                     highlight = 'prs07_enr',
                     ignore_offsets = True,
                     all_features = False,
                     pdfname = 'pp_features_enric_prs0.7.pdf')

l.fractions_barplot2(features = True,
                     highlight = 'wnz',
                     ignore_offsets = True,
                     all_features = False,
                     pdfname = 'pp_features_wrong_nonzero.pdf')

# ##

l.peak_ratio_score_threshold = 0.5
l.peak_ratio_score_bool()
l.fractions_barplot2(features = True,
                     highlight = 'prs1',
                     ignore_offsets = True,
                     all_features = False,
                     pdfname = 'pp_features_prs_new_0.5.pdf')

l.peak_ratio_score_best(best = 5)
l.fractions_barplot2(features = True,
                     highlight = 'prsb',
                     ignore_offsets = True,
                     all_features = False,
                     pdfname = 'pp_features_prs_new_best5.pdf')

l.fractions_barplot2(features = True,
                     highlight = 'slobb',
                     ignore_offsets = True,
                     all_features = False,
                     pdfname = 'pp_features_enric_wo_wrong.pdf')

l.logical_and(one = 'prs1', two = 'slobb', result = 'prs07_slo')

