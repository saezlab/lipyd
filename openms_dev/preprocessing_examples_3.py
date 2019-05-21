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
    resampleid = re.compile(r'\d+_\w+_([A-Z]\d{2})_(?:pos|neg).*')
    name = os.path.basename(name)
    return resampleid.match(name).groups()[0]


pp = msproc.PeakPickerHiRes(
    input_path = example_data_0,
    sample_id_method = sample_id_method,
)
pp.main()

ff = msproc.FeatureFindingMetabo(
    input_path = pp.output_path,
    sample_id = pp.sample_id,
)
ff.main()
