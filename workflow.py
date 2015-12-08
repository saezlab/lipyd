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

'''
Initializing from the beginning.
'''
#data, fnames, samples, csamples, pprofs = \
#    init_from_scratch(basedir, ltpdirs, pptablef, samplesf)

'''
Here we reload the workspace including all the data from the pickle:
'''
data, fnames, samples, csamples, samples_upper, pprofs = ltp.init_reinit(ltp.basedir)

'''
Setting the protein concentrations to zero in controls
out of protein profiles.
'''
ltp.zero_controls(samples_upper, pprofs)

'''
Reloading profiles to know the original values.
'''
pprofs_original = ltp.protein_profiles(ltp.ppsecdir, ltp.ppfracf)

'''
Plotting barcharts of protein profiles.
'''
ltp.fractions_barplot(samples_upper, pprofs, pprofs_original)

'''
Apply basic filters, and obtaining the valid features.
'''
ltp.apply_filters(data)
ltp.validity_filter(data)
valids = ltp.valid_features(data)

'''
Normalizing all the features.
'''
ltp.norm_all(valids)

'''
Correlation and similarity metrics between features and
protein concentration profile.
'''
ltp.profiles_corrs(valids, pprofs, samples_upper)

'''
Calculating ubiquity.
'''
ltp.sort_alll(valids, 'mz')
ltp.ubiquity_filter(valids)

'''
Lookup lipids for all the valid features.
'''
exacts, runtime = ltp.lipid_lookup_exact(valids, ltp.swisslipids_url)

'''
Looking up gold standard.
'''
stdpos = ltp.read_positives(ltp.basedir)
ltp.sort_alll(valids, 'mz')
for ltpname in stdpos.keys():
    ltp.true_positives(valids, stdpos, ltpname, 'pos')
    ltp.true_positives(valids, stdpos, ltpname, 'neg')

'''
Sensitivity, specificity and ROC curve
'''
score_perf = ltp.evaluate_scores(valids, stdpos.keys())
ltp.plot_score_performance(score_perf)
ltp.plot_roc(score_perf)

'''
Best hits according to combined Euclidean distance
and Goodman-Kruskal's gamma.
'''
ltp.best_gk_eu(valids)

'''
Exporting lists of 10 best features per LTP.
'''
ltp.best_table(valids, 'best10_positive.csv', 'pos')
ltp.best_table(valids, 'best10_negative.csv', 'neg')

pFragments = ltp.read_metabolite_lines('lipid_fragments_positive_mode.txt')
nFragments = ltp.read_metabolite_lines('lipid_fragments_negative_mode.txt')
ms2files = ltp.ms2_filenames(ltp.ltpdirs)
ms2map = ltp.ms2_map(ms2files)
ltp.ms2_main(valids, samples_upper, ms2map, pFragments, nFragments)

#################################################################
#########################################################
###############################
################
######
####
ltp.basic_filters(data, pprofs, samples, csamples)
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
ms2files = ltp.ms2_filenames(ltp.ltpdirs)
ms2map = ms2_map(ms2files)
ms2_main(ms2map, stage2_best, pFragments, nFragments)
ms2_main(ms2map, stage2_best_unknown, pFragments, nFragments)