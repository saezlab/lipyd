#!/usr/bin/env python
#-*- coding: utf-8 -*-

# Denes Turei
# turei.denes@gmail.com

import os
import imp
import re

import pyopenms as oms

import lipyd.msproc_alternative as msproc
import lipyd.settings as settings
import lipyd.common as common

example_data_dir = os.path.join(
    '/', 'home', 'denes', 'archive', 'ltp',
    'STARD10_invivo_raw', 'mzml', 'pos',
)
example_data_0 = os.path.join(
    example_data_dir,
    '150310_Popeye_MLH_AC_STARD10_A10_pos.mzML'
)

def reload():
    imp.reload(common)
    imp.reload(settings)
    imp.reload(msproc)


def sample_id_method(name):
    resampleid = re.compile(r'\d+_Popeye_MLH_AC_(\w+)_(?:pos|neg).*')
    name = os.path.basename(name)
    return resampleid.match(name).groups()[0]

def sample_filter_method(name):
    return 'MeOH' not in name and 'buffer_extracted' not in name

# single instance test: peak picking
pp = msproc.PeakPickerHiRes(
    input_path = example_data_0,
    sample_id_method = sample_id_method,
)
pp.main()

# single instance test: feature finding
ff = msproc.FeatureFindingMetabo(
    input_path = pp.output_path,
    sample_id = pp.sample_id,
)
ff.main()

# single instance test: map alignment
# ok, this does not make sense,
# the transformation should result identical map
ma = msproc.MapAlignmentAlgorithmPoseClustering(
    input_path = ff.output_path,
    reference_path = ff.output_path,
    sample_id = ff.sample_id,
)
ma.main()

# single instance test: MGF export
mgfe = msproc.MgfExport(
    input_path = pp.output_path,
    sample_id = pp.sample_id,
)
mgfe.main()

# pipeline test
wf = msproc.MSPreprocess(
    input_path = example_data_dir,
    stop = 'features_aligned',
    sample_id_method = sample_id_method,
    reference_sample = 'STARD10_A10',
)
# step by step
wf.peak_picking()
wf.feature_finding()
wf.map_alignment()
wf.feature_grouping()
wf.export()

# single call
wf = msproc.MSPreprocess(
    input_path = example_data_dir,
    stop = 'features_aligned',
    sample_id_method = sample_id_method,
    reference_sample = 'STARD10_A10',
)
wf.main()
