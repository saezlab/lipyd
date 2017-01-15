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

from future.utils import iteritems

#from common import *

#

from ltp import ltp2

#l = ltp2.LTP(recal_source = 'denes', only_marcos_fragments = False, marco_lipnames_from_db = False)
l = ltp2.Screening()


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

l.std_layout_tables_xlsx()
l.read_marco_standards()
l.standards_crosscheck(only_best = True)

### Diff between recent and old erroneous output

# l = ltp2.Screening(manual_ppratios_xls = 'Proteins_Overview_03b.xlsx', manual_ppratios_xls_cols = [2, 4, 8, 9])
# l.features_xls_diff('top_features_July2016', 'top_features')
l.features_xls_diff('top_features_2016Oct13_backup', 'top_features')
l.features_xls_diff('top_features_2016Oct14', 'top_features')


### PTPN9 unknown headgroups

ptpn9 = [338.1387, 784.6218, 756.5908, 518.3639, 730.5754, 758.6063, 520.3788, 492.3486]

from ltp import ltp2

ptpn9 = \
    [
        [338.143, 292.1328],           # ??
        [492.354, 255.2299],           # 16:0
        [518.370, 281.2452],           # 18:1
        [520.385, 283.2607],           # 18:0
        [730.584, 255.2298, 283.2610], # 16:0, 18:0
        [756.600, 281.2451, 255.2298], # 18:1, 16:0
        [756.600, 253.2142, 283.2604], # 16:1, 18:0
        [758.615, 283.2610, 255.2298], # 18:0, 16:0
        [784.622, 281.2452, 283.2607], # 18:1, 18:0
        [804.570, 281.2465, 255.2309]  # control: PC(18:1/16:0),
                                       # returns gly+P+choline (223.09734506418)
                                       # assuming formiate adduct
    ]

def get_headgroup_mass(m, ads = ['h', 'fo']):
    """
    Returns residual mass after calculating exact mass
    and removing fatty acid masses.
    
    :param list ads: Adducts assumed.
    """
    return \
        list(
            map(
            lambda ad:
                (
                    ad,
                    getattr(ltp2.Mz(m[0]), 'remove_%s' % ad)() - sum(m[1:]),
                ),
                ads
            )
        )

headgroups = \
    list(
        map(
            get_headgroup_mass,
            ptpn9
        )
    )

# resulted:
#[
    #[('h', 45.00292353322999),   ('fo', 1.0119971474307476)], # this does not make sense
    
    #[('h', 236.11682353323),     ('fo', 192.12589714743075)], # for 492.4, 518.4, 520.4:
    #[('h', 236.11752353322998),  ('fo', 192.12659714743074)], # headgroup has a mass
    #[('h', 236.11702353323),     ('fo', 192.12609714743076)], # either 236 or 192
    
    #[('h', 191.0859235332299),   ('fo', 147.09499714743072)], # 730.6 apparently does not
                                                               # fit in these series
    
    #[('h', 219.11782353323008),  ('fo', 175.1268971474309)],  # for 756.6, 758.6, 784.6:
    #[('h', 219.11812353323),     ('fo', 175.12719714743082)], # headgroup has a mass
    #[('h', 219.11692353322996),  ('fo', 175.12599714743078)], # either 219 or 175
    #[('h', 219.10882353322995),  ('fo', 175.11789714743077)],
    #[('h', 267.08532353323005),  ('fo', 223.09439714743087)]  # control, gly+P+choline mass correct
#]

def find_constitutions(
                    target,
                    elems = ['C', 'H', 'O', 'N', 'P', 'S'],
                    tolerance = 0.02,
                    maxcnts = None,
                    counts = None,
                    results = []
    ):
    c = int(np.ceil((target - 2) / ltp2.MolWeight('CH2').weight))
    h = ltp2.mas.mass['H']
    counts = dict(zip(elems, [0] * len(elems)))
    counts['C'] = c
    counts['H'] = int(np.round((target - ltp2.MolWeight('C%u' % c)) / h))
    
    
    maxcnts = \
        maxcnts or \
        list(map(lambda e: int(target / ltp2.MolWeight('%s' % e).weight), elems))
    counts = counts or [0] * (len(elems) - 1)
    const = ''.join(map(lambda e: '%s%u' % e, zip(elems, counts)))
    w = ltp2.MolWeight(const).weight
    if w >= target + tolerance
    for i in xrange(len(elems)):
        if counts[i] < maxcnts[i]:
            counts[i] += 1
            break
    

###

# output ppratio ranges:
    
old_ratios = copy.deepcopy(l.first_ratio)

with open('ppratio_ranges.csv', 'w') as f:
    for protein in sorted(l.first_ratio.keys()):
        r = l.first_ratio[protein]
        if len(r) > 1:
            sys.stdout.write('%s has more than 1 pairs of fractions\n' % protein)
        fr, rng = r.items()[0]
        f.write('%s\t%s:%s\t%.03f\t%.03f\n' % (protein, fr[0], fr[1], rng[0], rng[1]))

###

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
