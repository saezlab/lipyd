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
import time
import os
import sys
import copy
import cPickle as pickle

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
pprofs_original = copy.deepcopy(pprofs)
singles = ltp.one_sample(samples_upper)

lipnames = ltp.read_lipid_names(ltp.lipnamesf)
bindprop = ltp.read_binding_properties(ltp.bindpropf)

'''
Setting the protein concentrations to zero in controls
out of protein profiles.
'''
ltp.zero_controls(samples_upper, pprofs)

'''
Reloading profiles to know the original values.
'''
# pprofs_original = ltp.protein_profiles(ltp.ppsecdir, ltp.ppfracf, fnames)

'''
Apply basic filters, and obtaining the valid features.
'''
if data is not None:
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
Clustering features based on euclidean distance
'''
ltp.distance_matrix(valids, metrics = ['en'], with_pprof = True, 
    pprofs = pprofs, samples = samples_upper)

ltp.features_clustering(valids)
ltp.distance_corr(valids)

reload(ltp)
ltp.plot_heatmaps_dendrograms(valids, singles, pprofs, samples_upper,
    threshold = 0.05, threshold_type = 'percent', coloring = 'dist',
    ltps = ['STARD10'], fname = 'features_clustering_STARD10.pdf')

ltp.plot_heatmaps_dendrograms(valids, singles, ltps = ['STARD10'], 
    threshold = 0.05, threshold_type = 'percent', coloring = 'dist')


ltp.plot_heatmaps_dendrograms(valids, singles, ltps = ['STARD10'], 
    threshold = None, threshold_type = 'clsize', coloring = 'dist')


ltp.plot_heatmaps_dendrograms(valids, singles, ltps = ['STARD10'], 
    threshold = 30, threshold_type = 'incons', coloring = 'dist')

# run all with clsize
ltp.plot_heatmaps_dendrograms(valids, singles, pprofs, samples_upper,
    threshold = None, threshold_type = 'clsize', coloring = 'dist',
    save_selection = 'clcs10pct')

bs = ltp.fractions_barplot(samples_upper, pprofs, pprofs_original,
    features = True, valids = valids, 
    highlight = 'clcs10pct',
    highlight2 = False,
    all_features = False,
    pdfname = 'protein_profiles_features_cl_csize_10pct.pdf')

# run for all with 5% threshold
ltp.plot_heatmaps_dendrograms(valids, singles, pprofs, samples_upper, 
    threshold = 0.70, threshold_type = 'percent', coloring = 'dist',
    save_selection = 'cl70pct')

bs = ltp.fractions_barplot(samples_upper, pprofs, pprofs_original,
    features = True, valids = valids, 
    highlight = 'cl50pct',
    highlight2 = False,
    all_features = False,
    pdfname = 'protein_profiles_features_cl_maxdist_50pct.pdf')

# ##

'''
Calculating ubiquity. (variable name: `ubi`)
'''
ltp.sort_alll(valids, 'mz')
ltp.ubiquity_filter(valids)
timeit.timeit('ltp.ubiquity_filter(valids)', 
    setup = 'from __main__ import ltp, valids', number = 1)

'''
Lookup lipids for all the valid features. (variable name: `lip`)
'''
exacts = None
exacts, runtime = ltp.lipid_lookup_exact(valids, ltp.swisslipids_url, exacts = exacts, lipnames = lipnames)
ltp.negative_positive2(valids, lipnames)

ltp.ms1_headgroups(valids, lipnames, verbose = False)

ltp.headgroups_negative_positive(valids, 'ms1')

# the MS2 part
pFragments, pHgfrags, pHeadgroups = ltp.read_metabolite_lines('lipid_fragments_positive_mode.txt')
nFragments, nHgfrags, nHeadgroups = ltp.read_metabolite_lines('lipid_fragments_negative_mode.txt')
ms2files = ltp.ms2_filenames(ltp.ltpdirs)
ms2map = ltp.ms2_map(ms2files)
ltp.ms2_main(valids, samples_upper, ms2map, pFragments, nFragments)
ltp.ms2_headgroups(valids, pHgfrags, nHgfrags, pHeadgroups, nHeadgroups)
ltp.headgroups_negative_positive(valids, 'ms2')
#

ltp.feature_identity_table(valids)

ms1tab_coln, ms1tab = ltp.ms1_table(valids, lipnames)


ltp.ms1_table_html(valids, lipnames)
ltp.ms1_table_html_simple(valids, lipnames, include = 'cl5pct')
ltp.ms2_table_html_simple(valids, lipnames, include = 'cl5pct')

ltp.ms1_ms2_table_html_simple(valids, lipnames, include = 'cl70pct')

idlevels = {
    'All': ltp.identification_levels(valids, 'STARD10', 'PC'),
    'Best': ltp.identification_levels(valids, 'STARD10', 'PC', classif = 'cl50pct')
}

ltp.plot_identification_levels(idlevels, 'STARD10', 'PC')

idlevels = {
    'All': ltp.identification_levels(valids, 'LCN1', 'PC'),
    'Best': ltp.identification_levels(valids, 'LCN1', 'PC', classif = 'cl50pct')
}

ltp.plot_identification_levels(idlevels, 'LCN1', 'PC')

idlevels = {
    'All': ltp.identification_levels(valids, 'STARD2', 'PC'),
    'Best': ltp.identification_levels(valids, 'STARD2', 'PC', classif = 'cl50pct')
}

ltp.plot_identification_levels(idlevels, 'STARD2', 'PC')

## THIS IS NOT NEEDED ANY MORE ##

ltp.count_threshold_filter(valids, 'euv', threshold = 3.3, count = 10)
ltp.count_threshold_filter(valids, 'env', threshold = 0.01, threshold_type = 'best_fraction', count = 10)

ltp.count_threshold_filter(valids, 'env', threshold = 0.33, threshold_type = 'median_relative', count = 1000)

ltp.scores_plot(valids, score = 'env', singles = singles, pdfname = 'scores_median_33pct.pdf')

bs = ltp.fractions_barplot(samples_upper, pprofs, pprofs_original,
    features = True, valids = valids, 
    highlight = 'bool_env',
    highlight2 = False,
    all_features = True,
    pdfname = 'protein_profiles_features_all_median_33pct.pdf')

bs = ltp.fractions_barplot(samples_upper, pprofs, pprofs_original,
    features = True, valids = valids, 
    highlight = 'bool_env',
    highlight2 = False,
    all_features = False,
    pdfname = 'protein_profiles_features_median_33pct.pdf')

ltp.count_threshold_filter(valids, 'env', threshold = 0.33, threshold_type = 'mean_relative', count = 1000)

bs = ltp.fractions_barplot(samples_upper, pprofs, pprofs_original,
    features = True, valids = valids, 
    highlight = 'bool_env',
    highlight2 = False,
    all_features = True,
    pdfname = 'protein_profiles_features_all_mean_33pct.pdf')

bs = ltp.fractions_barplot(samples_upper, pprofs, pprofs_original,
    features = True, valids = valids, 
    highlight = 'bool_env',
    highlight2 = False,
    all_features = False,
    pdfname = 'protein_profiles_features_mean_33pct.pdf')

ltp.count_threshold_filter(valids, 'env', threshold = 0.050, threshold_type = 'fix', count = 100)

bs = ltp.fractions_barplot(samples_upper, pprofs, pprofs_original,
    features = True, valids = valids, 
    highlight = 'bool_env',
    highlight2 = False,
    all_features = True,
    pdfname = 'protein_profiles_features_all_fix_0.05.pdf')

bs = ltp.fractions_barplot(samples_upper, pprofs, pprofs_original,
    features = True, valids = valids, 
    highlight = 'bool_env',
    highlight2 = False,
    all_features = False,
    pdfname = 'protein_profiles_features_fix_0.05.pdf')

ltp.count_threshold_filter(valids, 'env', threshold = 0.10, threshold_type = 'fix', count = 100)

# relative to min
ltp.count_threshold_filter(valids, 'env', threshold = 10, threshold_type = 'relative', count = 100)

bs = ltp.fractions_barplot(samples_upper, pprofs, pprofs_original,
    features = True, valids = valids, 
    highlight = 'bool_env',
    highlight2 = False,
    all_features = True,
    pdfname = 'protein_profiles_features_all_fix_0.1.pdf')

bs = ltp.fractions_barplot(samples_upper, pprofs, pprofs_original,
    features = True, valids = valids, 
    highlight = 'bool_env',
    highlight2 = False,
    all_features = False,
    pdfname = 'protein_profiles_features_fix_0.1.pdf')

# ## ## ##

#pn_ratios = {}
#p_counts = {}
#n_counts = {}
#for l, d in valids.iteritems():
    #pn_ratios[l] = sum(d['pos']['bool_env']) / float(sum(d['neg']['bool_env']))
    #p_counts[l] = sum(d['pos']['bool_env'])
    #n_counts[l] = sum(d['neg']['bool_env'])
    #print l, sum(d['pos']['bool_env']), sum(d['neg']['bool_env'])

#for l, d in valids.iteritems():
    #for pn, tbl in d.iteritems():
        #tbl['env30'] = np.array([True] * 30 + [False] * (len(tbl['env']) - 20))


'''
Plotting barcharts of protein profiles.
'''
ltp.fractions_barplot(samples_upper, pprofs, pprofs_original)
ltp.fractions_barplot(samples_upper, pprofs, pprofs_original,
    features = False, valids = valids)


## THIS IS NOT NEEDED ANY MORE ##

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

## UNTIL HERE ##



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
