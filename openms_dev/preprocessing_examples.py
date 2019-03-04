#!/usr/bin/env python
#-*- coding: utf-8 -*-

# Denes Turei
# turei.denes@gmail.com

import imp

import pyopenms as oms

import lipyd.msproc as msproc

profile_mzml = (
    '/home/denes/archive/ltp/STARD10_invivo_raw/'
    'mzml/150310_Popeye_MLH_AC_STARD10_A10_pos.mzML'
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
