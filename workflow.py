#!/usr/bin/env python2
# -*- coding: utf-8 -*-

#
#  Copyright (c) 2014-2015 - EMBL-EBI
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
import os
import sys

from common import *
from ltp import *

#data, fnames, samples, csamples, pprops = \
#    init_from_scratch(basedir, ltpdirs, pptablef, samplesf)
data, fnames, samples, csamples, samples_upper, pprops = ltp.init_reinit(ltp.basedir)
ltp.apply_filters(data)
ltp.validity_filter(data)
valids = ltp.valid_features(data)
ltp.norm_all(valids)
ltp.profiles_corrs(valids, pprops, samples_upper)
ltp.sort_alll(valids, 'mz')
ltp.ubiquity_filter(valids)
exacts, runtime = ltp.lipid_lookup_exact(valids, ltp.swisslipids_url)
stdpos = ltp.read_positives(ltp.basedir)
ltp.sort_alll(valids, 'mz')
for ltpname in stdpos.keys():
    ltp.true_positives(valids, stdpos, ltpname, 'pos')
    ltp.true_positives(valids, stdpos, ltpname, 'neg')

score_perf = ltp.evaluate_scores(valids, stdpos.keys())
ltp.plot_score_performance(score_perf)
ltp.plot_roc(score_perf)
ltp.best_gk_eu(valids)
ltp.best_table(valids, 'best10_positive.csv', 'pos')
ltp.best_table(valids, 'best10_negative.csv', 'neg')
#################################################################
#########################################################
###############################
################
######
####
ltp.basic_filters(data, pprops, samples, csamples)
# stage0 :: feature filtering
stage0 = ltp.get_scored_hits(data)
# stage1 :: lipids
#pAdducts, nAdducts, lipids, runtime = lipid_lookup(stage0, pAdducts, nAdducts)
exacts, lipids, unknowns, runtime = ltp.lipid_lookup_exact(stage0, swisslipids_url)
# stage2 :: positive-negative
stage2 = negative_positive(lipids)
stage2_unknown = negative_positive(unknowns)
stage2_best = best_matches(lipids, stage2, minimum = 100000)
stage2_best_unknown = best_matches(unknowns, stage2_unknown, minimum = 100000)
stage2_best_all = best_matches((lipids, unknowns), (stage2, stage2_unknown), minimum = 100000)
# output
write_out(stage2_best, 'lipid_matches_nov25.csv')
write_out(np.vstack((stage2_best, stage2_unknown)), 'all_sorted_nov25.csv')
# evaluation
evaluate_results(stage0, stage2, lipids, samples_upper, 'f')
# LTP, num of lipid hits (min), num of positive-negative matching, 
# num of selected features [neg, pos], num of fractions in sample
[(l.upper(), p, len(j), len(stage2[l.upper()])) if j is not None else (l.upper(), p, 0, 0) \
    for l, i in lipids.iteritems() for p, j in i.iteritems()]
# processing MS2
pFragments = read_metabolite_lines('lipid_fragments_positive_mode.txt')
nFragments = read_metabolite_lines('lipid_fragments_negative_mode.txt')
ms2files = ms2_filenames(ltpdirs)
ms2map = ms2_map(ms2files)
ms2_main(ms2map, stage2_best, pFragments, nFragments)
ms2_main(ms2map, stage2_best_unknown, pFragments, nFragments)