#!/usr/bin/env python
#-*- coding: utf-8 -*-

# Denes Turei
# turei.denes@gmail.com

import os
import imp
import re

import pyopenms as oms

import lipyd.msproc as msproc

example_data_dir = os.path.join(
    '/', 'home', 'denes', 'archive', 'ltp',
    'STARD10_invivo_raw', 'mzml', 'pos',
)
example_data_0 = os.path.join(
    example_data_dir,
    '150310_Popeye_MLH_AC_STARD10_A10_pos.mzML'
)


def sample_id_method(name):
    
    resampleid = re.compile(r'\d+_\w+_([A-Z]\d{2})_(?:pos|neg)\.mzML')
    
    resampleid.match(name).groups()


msproc.PeakPickerHiRes(
    input_path = example_data_0,
    
)

#p = msproc.MSPreprocess(
    #profile_mzml = profile_mzml,
    #smooth_profile_data = True,
    #gaussian_smoothing_param = {
        #'gaussian_width': .7,
    #},
    #export_smoothed_profile = True,
    #console = False,
#)
#p.main()


p = msproc.MSPreprocess(
    profile_mzml = profile_mzml,
    smooth_profile_data = False,
    export_smoothed_profile = False,
    console = False,
)
p.main()


p = msproc.MSPreprocess(
    profile_mzml = profile_mzml,
    smooth_profile_data = True,
    gaussian_smoothing_param = {
        'gaussian_width': .1,
    },
    export_smoothed_profile = True,
    console = False,
)
p.main()
