#!/usr/bin/env python2
# -*- coding: utf-8 -*-

#
#  Copyright (c) 2015-2016 - EMBL-EBI
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

#

from ltp import ltp2

l = ltp2.LTP()

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
l.recalibration()

'''
Lookup lipids for all the valid features. (variable name: `lip`)
'''
l.ms1()
l.ms2()

l.std_layout_tables_xls()

recal_marco = l.ltps_drifts
l.recal_source = 'denes'
l.recalibration()

with open('recalibration.txt', 'w') as f:
    for mode in ['neg', 'pos']:
        for protein in sorted(l.ltps_drifts.keys()):
            for frac in sorted(l.ltps_drifts[protein][mode].keys()):
                if frac not in recal_marco[protein.upper()][mode]:
                    print protein, mode, frac
                else:
                    f.write(
                        '%s\t%s\t%s\t%.03f\t%.03f\n' % (
                            protein.upper(), mode, frac,
                            recal_marco[protein.upper()][mode][frac],
                            l.ltps_drifts[protein][mode][frac]
                        )
                    )

l.identify()

'''
Correlation and similarity metrics between features and
protein concentration profile.
'''
l.profiles_corrs()
l.profiles_corrs(pprofs = '15')
l.profiles_corrs(pprofs = '45')

l.ms2_scans_identify()
l.find_known_binders()
l.known_binders_as_standard()

l.protein_peak_ratios()
l.intensity_peak_ratios()
l.peak_ratio_score()
l.peak_ratio_score_bool()

l.evaluate_scores()


'''
Clustering features based on euclidean distance
'''
l.distance_matrix(metrics = ['en'], with_pprof = True)
l.features_clustering()
l.plot_heatmaps_dendrograms(
    threshold = 0.70, threshold_type = 'percent', coloring = 'dist',
    save_selection = 'cl70pct')

'''
Calculating ubiquity. (variable name: `ubi`)
'''
l.sort_alll('mz')
l.ubiquity_filter()

'''
Exporting HTML tables
'''
l.ms1_table_html()
l.ms1_table_html_simple()
l.ms2_table_html_simple()
l.ms1_ms2_table_html_simple()
l.features_table()


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


idlevels = {
    'All': ltp.identification_levels(valids, 'STARD10', 'PC'),
    'Best': ltp.identification_levels(valids, 'STARD10', 'PC',
        classif = 'cl50pct')
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


####
ms2result = {
    'STARD10': {
        'neg': ['PC', 'PE'],
        'pos': ['PC']
    },
    'BPI': {
        'neg': ['PC', 'PI'],
        'pos': ['PE']
    },
    'LCN1': {
        'neg': ['PC'],
        'pos': []
    },
    'BPIFB2': {
        'neg': ['PC'],
        'pos': []
    },
    'GM2A': {
        'neg': ['PC', 'PI', 'PE', 'PS', 'PG'],
        'pos': ['PC', 'PE']
    },
    'STARD2': {
        'neg': ['PC'],
        'pos': []
    },
    'SEC14L2': {
        'neg': ['PI'],
        'pos': []
    },
    'ORP9': {
        'neg': ['PS'],
        'pos': ['Cer']
    },
    'GLTPD1': {
        'neg': ['SM'],
        'pos': []
    },
    'RBP1': {
        'neg': ['PG'],
        'pos': []
    },
    'SEC14L6': {
        'neg': [],
        'pos': ['PE', 'SM']
    },
    'STARD11': {
        'neg': [],
        'pos': ['Cer']
    },
    'SEC14L1': {
        'neg': [],
        'pos': ['Cer']
    },
    'FABP4': {
        'neg': [],
        'pos': ['Cer']
    },
    'GLTP': {
        'neg': [],
        'pos': ['Cer']
    }
}

ms2i = l.ms2identities_summary()

for protein, d in ms2i.iteritems():
    for mode, lst in d.iteritems():
        ms2i[protein][mode] = set(ms2i[protein][mode])

for protein, d in ms2result.iteritems():
    for mode, lst in d.iteritems():
        ms2i[protein][mode] = ms2i[protein][mode] | set(lst)
