#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `lipyd` python module
#
#  Copyright (c) 2015-2019 - EMBL
#
#  File author(s): Dénes Türei (turei.denes@gmail.com)
#
#  Distributed under the GNU GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://denes.omnipathdb.org/
#

# from pyteomics import mzml

try:
    import pyopenms as oms
except:
    pass

mzml_example = ('/home/denes/archive/ltp/STARD10_invivo_raw/mzml/'
    '150310_Popeye_MLH_AC_STARD10_A10_pos.mzML')

class Reader(object):
    """
    Reads and preprocesses mzML files using the OpenMS library.
    """
    
    def __init__(self, fname):
        """
        Parameters
        ----------
        fname : str
            File name of the input mzML file.
        """
        
        self.fname = fname
        
    def mgf(self):
        """
        Opens an mzML file.
        """
        
        options = oms.PeakFileOptions()
        options.setMSLevels([2])
        self.mzml = oms.MzMLFile()
        self.mzml.setOptions(options)
        self.exp = oms.MSExperiment()
        self.mzml.load(self.fname, self.exp)
        self.feature_finder = oms.FeatureFinder()
        self.ffname = 'centroided'
        self.features = oms.FeatureMap()
        self.seeds = oms.FeatureMap()
        self.params = oms.FeatureFinder().getParameters(self.ffname)
        self.feature_finder.run(
            self.ffname,
            self.exp,
            self.features,
            self.params,
            self.seeds,
        )
        self.features.setUniqueIds()
