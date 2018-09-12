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

from future.utils import iteritems

import os
import collections

from lipyd import common

_defaults = {
    # The absolute root directory.
    # This should not be necessary, why is it here?
    'path_root': '/',
    # The basedir for every files and directories in the followings.
    'basedir': os.getcwd(),
    # If None will be the same as ``basedir``.
    'data_basedir': None,
    # the overview table of Enric's screening
    # we read the protein containing and the
    # highest fractions from here
    'invitro_fractionsf': 'LTPsceenprogres_v07d.xlsx',
    'invivo_fractionsf': 'control_sample.csv',
    # list of potentially problematic fractions
    'wrongfracsf': 'wrong_fractions.csv',
    # File with offsets of fractions from SEC in ml.
    # This we used at Antonella, but at Enric, these values
    # are different for each protein, and contained by the
    # individual SEC profile files.
    'ppfracf': 'fractions.csv',
    # Directory with the SEC absorbance profiles for each protein.
    'ppsecdir': 'SEC_profiles',
    # The first section in the SEC profile filenames is the protein
    # name or the second:
    'sec_filenames_protein_first': True,
    # This is a typo in some of the PEAK table headers.
    # If this cause trouble just turn it off.
    'fix_fraction_name_ab_typo': True,
    # Directory with the SDS PAGE protein quantities for some
    # proteins.
    'gelprofdir': 'gel_profiles',
    # at Antonella, these proteins come before our fractions, or
    # maybe in A9, so the further fractions we can consider void
    # and we can use their average to reduce the background of the
    # absorbances of the same fractions of other proteins.
    'background_proteins': set(['BNIPL', 'OSBP', 'SEC14L1']),
    # constant fraction layout in Antonella's screening
    # not used at Enric, is kind of "deprecated"
    'fixed_fraction_layout': False,
    'fracs': ['A9', 'A10', 'A11', 'A12', 'B1'],
    'all_fracs': ['A5', 'A6', 'A7', 'A8', 'A9',
                    'A10', 'A11', 'A12', 'B1'],
    # constant fraction layout in Antonella's screening
    # not used at Enric, is kind of "deprecated"
    # uppercase version with leading zeros
    'fracsU': ['A09', 'A10', 'A11', 'A12', 'B01'],
    'fr_offsets': {
        'invitro': [0.0],          # at Enric this is [0.0],
        'invivo':  [0.010, 0.045], # at Antonella [0.010, 0.045]
    }, 
    'abs_cols': [1], # this can be [1, 3, 5] if we want to read
                        # also the 260 and 215 nm UV profiles...
    'fraclims_from_sec': True, # read the fraction limits from
                                # the same files as SEC profiles
                                # or simply from `fractions.csv`
    'pp_do_correction': False, # bypass corrections
    # determine protein containing fractions from the absorbances
    # or read from separate file
    'pcont_fracs_from_abs': {
        'invivo': False,
        'invitro': True,
    },
    # table with manually set protein ratios among
    # many other columns
    'manual_ppratios_xls': 'Proteins_Overview_05.xlsx',
    # Columns to read from the manual ppratios XLS file.
    # These were different at Antonella and Enric, so we
    # need to set them here.
    #'manual_ppratios_xls_cols': [2, 4, 8, 9], # 03b
    'manual_ppratios_xls_cols': [3, 6, 10, 11], # 05
    # Read externally determined protein peak ratios from file
    # these were provided by Marco and used for Antonella`s
    # data analysis
    'use_manual_ppratios': {
        'invitro': False,
        'invivo': True,
    },
    'protein_containing_fractions_invivo':
        'protein_containing_fractions_invivo.tsv',
    'protein_containing_fractions_invivo':
        'protein_containing_fractions_invitro.tsv',
}

in_basedir = [
    'invivo_fractionsf', 'invitro_fractionsf', 'ppfracf', 'seqfile',
    'pptablef', 'lipnamesf', 'bindpropf', 'metabsf',
    'featurescache',
    'auxcache', 'stdcachefile', 'validscache', 'marco_dir',
    'abscache', 'pptable_file', 'recalfile', 'manual_ppratios_xls',
    'manualdir', 'ltplistf', 'flimcache', 'ppsecdir', 'gelprofdir']

in_datadir = {
    'pfragmentsfile', 'nfragmentsfile', 'lipnamesf', 'mgf_example',
    
}

def reset_all():
    
    settings = collections.namedtuple('Settings', list(_defaults.keys()))
    
    for k in _defaults.keys():
        
        val = getattr(defaults, k)
        
        if k in in_basedir:
            val = os.path.join(_defaults['basedir'], val)
        
        if k in in_datadir:
            val = os.path.join(_defaults['datadir'], val)
        
        setattr(settings, k, val)
    
    globals()['settings'] = settings

def setup(**kwargs):
    
    for param, value in iteritems(kwargs):
        
        setattr(settings, param, value)

def get(param):
    
    if hasattr(settings, param):
        
        return getattr(settings, param)

def get_default(param):
    
    if hasattr(defaults, param):
        
        return getattr(defaults, param)

def reset(param):
    
    setup(param, get_default(param))

defaults = common._const()

for k, v in iteritems(_defaults):
    
    setattr(defaults, k, v)

reset_all()
