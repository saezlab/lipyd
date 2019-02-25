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

p = msproc.MSPreprocess(
    profile_mzml = profile_mzml,
    smooth_profile_data = True,
)
p.main()
