#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `lipyd` python module
#
#  Copyright (c) 2015-2017 - EMBL
#
#  File author(s):
#  Dénes Türei (turei.denes@gmail.com)
#  Igor Bulanov
#
#  Distributed under the GNU GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://denes.omnipathdb.org/
#

from future.utils import iteritems

import os
import copy
import collections

import lipyd.common as common
import lipyd.lipproc as lipproc

_defaults = {
    # The absolute root directory.
    # This should not be necessary, why is it here?
    'path_root': '/',
    # The basedir for every files and directories in the followings.
    'basedir': os.getcwd(),
    # If None will be the same as ``basedir``.
    'data_basedir': None,
    # List of the shared directories containing most of the data.
    # Due to their large size these are kept on separate partition.
    # Lists are joined by `os.path.join`.
    # E.g. two shared folders can be defined as
    # `[['my', 'huge', 'partition'], [`another`, `huge`, `disk`]]`.
    'datadirs': [['share']],
    # the overview table of Enric's screening
    # we read the protein containing and the
    # highest fractions from here
    'fractionsf': 'LTPsceenprogres_v07d.xlsx',
    # list of potentially problematic fractions
    'wrongfracsf': 'wrong_fractions.csv',
    # File with offsets of fractions from SEC in ml.
    # This we used at Antonella, but at Enric, these values
    # are different for each protein, and contained by the
    # individual SEC profile files.
    'ppfracf': 'fractions.csv',
    # Directory with the SEC absorbance profiles for each protein.
    'ppsecdir': 'SEC_profiles',
    # background in SEC
    'sec_background': 'ARHGAP1',
    # The first section in the SEC profile filenames is the protein
    # name or the second:
    'sec_filenames_protein_first': True,
    # This is a typo in some of the PEAK table headers.
    # If this cause trouble just turn it off.
    'fix_fraction_name_ab_typo': True,
    # Directory with the SDS PAGE protein quantities for some
    # proteins.
    'gelprofdir': 'gel_profiles',
    # The directory containing the standards in mzML format.
    # These are used for calculation of recalibration.
    'stddir': 'Standards_mzML format',
    # Directory with manually processed files. These are originally
    # output of this module, then have been manually processed, and
    # can be read again and compared/further analysed.
    'manualdir': 'Processed_files',
    # ending of processed files
    'manualend': 'final.xlsx',
    # For recalibration we read the dates of runs from here.
    'seqfile': 'Sequence_list_LTP_screen_2015.csv',
    # This is an output file to export the absorbance based
    # protein quantities for all proteins and fractions.
    'pptablef': 'proteins_by_fraction.csv',
    # Defines abbreviations of each lipid names and the keywords
    # to identify these in SwissLipids and LipidMaps.
    'lipnamesf': 'lipid_names_v2.csv',
    # Literature curated data about known binding properties of LTPs.
    'bindpropf': 'binding_properties.csv',
    # Lipid classes properties and database IDsb
    'lipipropf': 'lipid_properties.csv',
    # The file with recalibration values from Marco.
    'recalfile': 'Recalibration_values_LTP.csv',
    # If the recalibration performed by this module, we read the
    # expected values of the standards from this file.
    'metabsf': 'Metabolites.xlsx',
    # Simple list with lipid binding domains, names and UniProt
    # IDs of all lipid transfer proteins.
    'ltplistf': 'ltplist.csv',
    # Tolerate numpy warnings, or raise them as errors.
    # Useful for debugging.
    'tolerate_numpy_warnings': True,
    # URLs of external databases.
    # SwissLipids: calculated masses of hundreds of thousands
    # of lipid species.
    'swisslipids_url': 'http://www.swisslipids.org/php/'\
        'export.php?action=get&file=lipids.csv',
    # Experimentally verified masses of few tens of thousands lipids.
    'lipidmaps_url': 'http://www.lipidmaps.org/resources/downloads/'\
        'LMSDFDownload12Dec17.tar.gz',
    # The filename to use after extracting the archice above.
    'lipidmaps_fname': 'LMSDFDownload12Dec17/'\
        'LMSDFDownload12Dec17FinalAll.sdf',
    # ComPPI: protein subcellular localization data.
    'comppi_url': 'http://comppi.linkgroup.hu/downloads',
    # Gene Ontology.
    'goa_url': 'ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/%s/'\
        'goa_%s.gaf.gz',
    # Gene Ontology.
    'quickgo_url': 'http://www.ebi.ac.uk/QuickGO/GAnnotation?format=tsv&'\
        'limit=-1%s&termUse=%s&tax=%u&col=proteinID,goID,goName,aspect',
    # PubChem webservice URL to query molecule properties
    # see details here:
    # https://pubchem.ncbi.nlm.nih.gov/help.html#Glossary
    'pubchem_url': ('https://pubchem.ncbi.nlm.nih.gov/'
        'rest/pug/compound/cid/%s/property/'
        'TPSA,' # polar surface area in square angstroms
        'XLogP,' # water octanol partitioning coefficient
        'Complexity,' # Bertz/Hendrickson/Ihlenfeldt formula
        'HBondDonorCount,' # O, N, P, S with hydrogene
        'HBondAcceptorCount,' # O, N, P, S with lone pair
        'HeavyAtomCount,' # non-H atom count
        'RotatableBondCount/XML'), # rotatable bondsgg
    # Manually curated data by Charlotte about localization
    # and membrane composition.
    'localizationf': 'subcellular_localisation_and_binding.xlsx',
    'membranesf': 'membranes_lipid_composition.xlsx',
    # Original MS2 fragment lists manually compiled by Marco.
    # One for positive and one for negative mode.
    'pfragmentsfile': 'lipid_fragments_positive_mode_v11d.txt',
    'nfragmentsfile': 'lipid_fragments_negative_mode_v11d.txt',
    # This is a file to output protein profiles to.
    # Serves only the purpose of checking for errors.
    'pptable_file': 'protein_profiles.txt',
    # Cache files. At certain points this module saves and reloads
    # data from pickles, so it don't need to reprocess everything
    # every time. Although if something changed, or something
    # goes wrong, you might need to delete these to ensure
    # everything processed by the most recent code.
    'featurescache': 'features.pickle',
    'pprofcache': 'pprofiles_raw.pickle',
    'abscache': 'absorbances.pickle',
    'flimcache': 'fraclims.pickle',
    # The directory with all MS2 MGF files. If not set, the MGF files
    # will be searched under the directory of each protein.
    'ms2dir': 'MGFfiles',
    # Directory with manually processed `golden standards`
    # from Marco.
    'marco_dir': 'marco',
    # table with manually set protein ratios among
    # many other columns
    'manual_ppratios_xls': 'Proteins_Overview_05.xlsx',
    # Columns to read from the manual ppratios XLS file.
    # These were different at Antonella and Enric, so we
    # need to set them here.
    #'manual_ppratios_xls_cols': [2, 4, 8, 9], # 03b
    'manual_ppratios_xls_cols': [3, 6, 10, 11], # 05
    'auxcache': 'save.pickle',
    'stdcachefile': 'calibrations.pickle',
    'validscache': 'valids.pickle',
    # Write a log about MS2 processing to this file.
    'ms2log': 'ms2identities.log',
    # use the items only of these levels from SwissLipids.
    # Useful to avoid flooding data with redundant subspecies.
    'swl_levels': ['Species'],
    # the source of recalibration: either from Marco, or calculated
    # by this software by processing raw data from the standards.
    'recal_source': 'marco',
    # at Antonella, these proteins come before our fractions, or
    # maybe in A9, so the further fractions we can consider void
    # and we can use their average to reduce the background of the
    # absorbances of the same fractions of other proteins.
    'background_proteins': set(['BNIPL', 'OSBP', 'SEC14L1']),
    # the minimum total intensities (average areas) of MS1 features
    # in negative and positive mode. Below these values features are
    # not highlighted and not included on the `best` sheet in the
    # output tables unless they have MS2 data.
    'aa_threshold': {
        'neg': 30000.0,
        'pos': 150000.0
    },
    'fr_offsets': [0.0],          # at Enric this is [0.0],
                                    # at Antonella [0.010, 0.045]
    'abs_cols': [1], # this can be [1, 3, 5] if we want to read
                        # also the 260 and 215 nm UV profiles...
    'fraclims_from_sec': True, # read the fraction limits from
                                # the same files as SEC profiles
                                # or simply from `fractions.csv`
    'pp_do_correction': False, # bypass corrections
    'pcont_fracs_from_abs': False, # determine protein containing
                                    # fractions from the absorbances
                                    # or read from separate file
    # allow features to have missing values in first or last
    # protein containing fractions
    'permit_profile_end_nan': True,
    # This is the maximum extent we can broaden the band
    # when selecting the features with intensity ratios
    # closest to the protein ratio. This is used as a
    # ratio-like limit, e.g. if the number here is 0.25,
    # and the protein ratio is 1.0, then the lower limit
    # will be 1.0 * 0.25 = 0.25 and the upper
    # 1.0 / 0.25 * 1.0 = 4.0.
    'peak_ratio_score_max_bandwidth': 0.25,
    'peak_ratio_score_optimal_bandwidth': 0.5,
    # This is the preferred number of features to
    # calculate the mean and SD in protein ratio score.
    # It means we try to select this number of features
    # with intensity ratio closest to the protein ratio,
    # and use their mean and SD to estimate the fit of 
    # other features. The lower the number the closer
    # features will be selected, but it should be enough
    # large to make the SD a meaningful metric. Similarly
    # we can set a minumum population, below this it
    # really does not make sense to calculate an SD.
    'peak_ratio_score_optimal_population': 10,
    'peak_ratio_score_still_good_population': 7,
    'peak_ratio_score_minimal_population': 4,
    # Use the adaptive method at the peak ratio score
    # calculation. This means to iteratively broaden the
    # band at selecting the features with intensity ratios
    # closest to the protein ratio, either until we have
    # the optimal population of features, or we reach the
    # maximum bandwith. If no features found this way,
    # a warning message will be displayed.
    'adaptive_peak_ratio_score': False,
    # Read externally determined protein peak ratios from file
    # these were provided by Marco and used for Antonella`s
    # data analysis
    'use_manual_ppratios': False,
    # omg, I don't remember what it is
    'use_last_ratio': False,
    # for calculation of peak ratio score, take all intensity ratios
    # within this range of tolerance around the protein ratio.
    # E.g if the protein ratio is 0.6, then the lower threshold
    # will be 0.6 * 0.5 = 0.3, and the upper 1 / 0.5 * 0.6 = 1.2;
    # Then the standard deviation of all intensity ratios within
    # this range is calculated from both positive and negative
    # features, and for each feature and for each pairs of fractions
    # the difference between its intensity ratio divided by the SD
    # results the `peak ratio score`. This is calculated for all
    # pairs of protein containing fractions, and we take their
    # mean if more than 2 fractions are available.
    'peak_ratio_range': 0.25,
    # The threshold below the peak ratio scores considered good.
    # We set a cut-off e.g. for highlighting with green in the
    # output tables, but otherwise it is a continuous measure
    # of goodness.
    'peak_ratio_score_threshold': 1.0,
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
    'pp_minratio': 3,
    # at Antonella's screening a fraction which is considered
    # void, so adjusting the zero of the absorbance profiles
    # to this fraction
    # see method ``pp_baseline_correction()``
    'basefrac': 'A5',
    # the tolerance when looking up MS1 m/z values in databases
    # ppm
    'ms1_tolerance': 10,
    # the tolerance when looking up MS2 fragment masses ppm
    'ms2_tolerance': 100,
    # tolarance at matching MS1 peaks against precursors in
    # mgf files with MS2 spectra
    'precursor_match_tolerance': 50,
    # the tolerance at identifying features in standards ppm
    'std_tolerance': 20,
    # MS2 precursors must have their charges determined
    'ms2_precursor_charge': None,
    # Use only fragments from Marco's lists (note: the series)
    # are still programmatically generated and their m/z values
    # calculated automatically), or use an extended fragment list
    # constructed by Denes (note: this might contain more false
    # positives, i.e. fragments which do not occure, or are iso-
    # baric with others, resulting unnecessary noise; but some-
    # times these annotations still might help to have a clue
    # what is there in cases where otherwise you would see only
    # unknown fragments)
    'only_marcos_fragments': True,
    # Expect only certain adducts at certain lipid categories
    # or assume any lipid may form any type of adduct.
    'adducts_constraints': False,
    # In output tables construct the lipid names programmatically
    # or copy the ones from the database.
    'marco_lipnames_from_db': True,
    # Use the average area values from the `PEAK` software
    # output, or recalculate them here using the intensity
    # values of features across all fractions
    'use_original_average_area': True,
    # Overwrite UV absorbance based protein quantities with
    # those manually acquired from SDS PAGE when these are
    # available
    'use_gel_profiles': True,
    # Use slope_filter (by Enric) for the final
    # selection of profiles
    'slope_profile_selection': True,
    # In slope_filter, if there are 2 highest fractions
    # with approximately equal protein content,
    # we use this range of tolerance to filter
    # the intensity ratios between them. E.g if this number is
    # 0.75, ratios between 0.75 and 1.33 will be accepted.
    'slope_equal_fractions_ratio': 0.75,
    # When we have only the descending parts of the profiles
    # in the measured fractions, we set this additional
    # constraint to avoid include everything descending
    # and select those which still fit better the protein
    'slope_descending_slope_diff_tolerance': 0.8,
    # Do not use the fractions marked as wrong
    # at comparing protein and intensity profiles
    'filter_wrong_fractions': True,
    # The MS2 retention time values must be within the RT
    # range of the feature detected in MS1, otherwise
    # will be dropped
    'ms2_rt_within_range': False,
    # Consider the MS2 scans from only those fractions
    # containing the protein, or from all available fractions
    'ms2_only_protein_fractions' : False,
    # Above this threshold we consider the MS2 spectrum to not
    # belong to the protein and highlight with red in the
    # output tables.
    'deltart_threshold': 0.5,
    # Don't know what it is for
    'uniprots': None,
    # example files
    'mgf_neg_examples': 'neg_examples.mgf',
    'mgf_pos_examples': 'pos_examples.mgf',
    'peaks_example': 'peaks_example.csv',
    'peaks_gltpd1_invitro': 'peaks_example_gltpd1_invitro_50.csv',
    'peaks_gltpd1_invivo': 'peaks_example_gltpd1_invivo_50.csv',
    'sec_gltpd1_invitro': 'SEC_GLTPD1_invitro.asc',
    'sec_gltpd1_invivo': 'SEC_GLTPD1_invivo.xls',
    'sec_xls_example': 'SEC_xls_example.xls',
    'sec_unicorn_example': 'SEC_asc_unicorn_example.asc',
    # an MFQL example file
    'mfql_example': 'Neg_bovine_heart_PE.mfql',
    # logarithm base for matching chain fragment intensities
    # this determines the tolerance and represents fold difference
    'chain_fragment_instensity_ratios_logbase': 1.5,
    # if no expected ratios provided still require the ratios to be even
    # at glycero(phospho)lipids?
    'even_chain_fragment_intensity_ratios_gl_gpl': True,
    # if no expected ratios provided still require the ratios to be even
    # at aphingolipids?
    'even_chain_fragment_intensity_ratios_sl': False,
    # at MS2 identification, add chain details
    # (fragment rank, intensity and type) to the MS2Identity object
    # turning this off makes MS2 spectra analysis faster
    'ms2_scan_chain_details': True,
    # Method names to convert between adduct and exact masses
    'ad2ex': {
        1: {
            'pos': {
                '[M+H]+': 'remove_h',
                '[M+NH4]+': 'remove_nh4',
                '[M+Na]+': 'remove_na',
                '[M-H2O+H]+': 'add_oh',
            },
            'neg': {
                '[M-H]-': 'add_h',
                '[M+HCOO]-': 'remove_fo',
            }
        },
        2: {
            'pos': {},
            'neg': {
                '[M-2H]2-': 'add_2h'
            }
        },
        3: {
            'pos': {},
            'neg': {
                '[M-3H]3-': 'add_3h'
            }
        }
    },
    # method names to convert between exact and adduct masses
    'ex2ad': {
        1: {
            'pos': {
                '[M+H]+': 'add_h',
                '[M+NH4]+': 'add_nh4',
                '[M+Na]+': 'add_na',
                '[M-H2O+H]+': 'remove_oh',
            },
            'neg': {
                '[M-H]-': 'remove_h',
                '[M+HCOO]-': 'add_fo',
            }
        },
        2: {
            'pos': {},
            'neg': {
                '[M-2H]2-': 'remove_2h'
            }
        },
        3: {
            'pos': {},
            'neg': {
                '[M-3H]3-': 'remove_3h'
            }
        }
    },
    # metrics to use for determining similarity of protein
    # and intensity profiles; these were attempts, and not
    # used any more, are "deprecated"
    'metrics': [
        ('Kendall\'s tau', 'ktv', False),
        ('Spearman corr.', 'spv', False),
        ('Pearson corr.', 'pev', False),
        ('Euclidean dist.', 'euv', True),
        ('Robust corr.', 'rcv', False),
        ('Goodman-Kruskal\'s gamma', 'gkv', False),
        ('Difference', 'dfv', True)
    ],
    # adducts used by default
    'adducts_default': {
        'neg': {
            1: {
                '[M-H]-',
                '[M+HCOO]-',
            },
        },
        'pos': {
            1: {
                '[M+H]+',
                '[M+NH4]+',
                '[M+Na]+',
            },
        }
    },
    # additional constraints for adduct lookups at various species
    # e.g. by default `[M-H2O+H]+` is not used but at Vitamin A we use it:
    'adduct_constraints': {
        'pos': {
            lipproc.Headgroup(main = 'VA'): {
                '[M+H]+',
                '[M+NH4]+',
                '[M+Na]+',
                '[M-H2O+H]+',
            },
            lipproc.Headgroup(main = 'Cer', sub = ('Hex',)): {
                '[M+H]+',
                '[M+NH4]+',
                '[M+Na]+',
                '[M-H2O+H]+',
            },
            lipproc.Headgroup(main = 'Cer', sub = ('Hex2',)): {
                '[M+H]+',
                '[M+NH4]+',
                '[M+Na]+',
                '[M-H2O+H]+',
            },
            lipproc.Headgroup(main = 'Cer', sub = ('SHex',)): {
                '[M+H]+',
                '[M+NH4]+',
                '[M+Na]+',
                '[M-H2O+H]+',
            },
            lipproc.Headgroup(main = 'Cer', sub = ('SHex2',)): {
                '[M+H]+',
                '[M+NH4]+',
                '[M+Na]+',
                '[M-H2O+H]+',
            },
        },
        'neg': {
        }
    },
    'font_family': ['Helvetica Neue LT Std', 'sans'],
    'font_variant': None,
    'font_style': None,
    'font_stretch': None,
    'font_size': 11,
    'font_sizes': {
        'axis_label': 1.,
        'ticklabel': .8,
        'legend_title': 1.,
        'legend_label': .8,
        'title': 1.2,
        'annotation': .8,
    },
    'figsize': (3, 4),
    'spectrum_plot_figsize': (9, 5),
    'spectrum_plot_xlab': 'm/z',
    'cachedir': None,
    # use only MS2 scans within the RT range of the feature
    'ms2_check_rt': True,
    'ms_preproc_wd': 'lipyd_ms_preproc',
}

in_basedir = [
    'fractionsf',
    'ppfracf',
    'seqfile',
    'pptablef',
    'lipnamesf',
    'bindpropf',
    'metabsf',
    'featurescache',
    'auxcache',
    'stdcachefile',
    'validscache',
    'marco_dir',
    'abscache',
    'pptable_file',
    'recalfile',
    'manual_ppratios_xls',
    'manualdir',
    'ltplistf',
    'flimcache',
    'ppsecdir',
    'gelprofdir',
    'ms_preproc_wd',
]

in_datadir = {
    'pfragmentsfile',
    'nfragmentsfile',
    'lipnamesf',
    'mgf_example',
    'peaks_example',
    'mfql_example',
    'sec_xls_example',
    'sec_unicorn_example',
    'peaks_gltpd1_invitro',
    'peaks_gltpd1_invivo',
    'sec_gltpd1_invitro',
    'sec_gltpd1_invivo',
}

in_mgfdir = {
    'mgf_neg_examples', 'mgf_pos_examples',
}


def reset_all():
    """ """
    
    settings = collections.namedtuple('Settings', list(_defaults.keys()))
    
    for k in _defaults.keys():
        
        val = getattr(defaults, k)
        
        if k in in_datadir:
            val = os.path.join(common.ROOT, 'data', val)
        
        if k in in_mgfdir:
            val = os.path.join(common.ROOT, 'data', 'ms2_examples', val)
        
        setattr(settings, k, val)
    
    if settings.cachedir is None:
        
        settings.cachedir = os.path.join(
            os.path.expanduser('~'),
            '.lipyd',
            'cache',
        )
    
    globals()['settings'] = settings


def setup(**kwargs):
    """

    Parameters
    ----------
    **kwargs :
        

    Returns
    -------

    """
    
    for param, value in iteritems(kwargs):
        
        setattr(settings, param, value)

def get(param):
    """

    Parameters
    ----------
    param :
        

    Returns
    -------

    """
    
    if hasattr(settings, param):
        
        return getattr(settings, param)

def get_default(param):
    """

    Parameters
    ----------
    param :
        

    Returns
    -------

    """
    
    if hasattr(defaults, param):
        
        return getattr(defaults, param)

def reset(param):
    """

    Parameters
    ----------
    param :
        

    Returns
    -------

    """
    
    setup(param, get_default(param))

defaults = common._const()

for k, v in iteritems(_defaults):
    
    setattr(defaults, k, v)

reset_all()
