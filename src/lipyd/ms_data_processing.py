#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `lipyd` python module
#
#  Copyright (c) 2015-2019 - EMBL
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

import pyopenms as pms


class MSPreprocess(object):
    """
    Main class for ms data processing using pyopenms library.

    Parameters
    ----------
    
    """
    
    
    def __init__(
            self,
            raw_mzml = None,
            peaks_file = None,
            features_file = None,
            seeds = None,
            fh = None,
            options = None,
            picked_input_map = None,
            ff = None,
            features = None,
            name = None,
            params = None,
            raw_map = None
        ):
        
        self.features_file = features_file
        self.raw_mzml = raw_mzml
        self.peaks_file = peaks_file
        self.raw_map = pms.MSExperiment()
        self.picked_out_map = pms.MSExperiment()
        self.seeds = pms.FeatureMap()
        self.fh = pms.MzMLFile()
        self.options = pms.PeakFileOptions()
        self.picked_input_map = pms.MSExperiment()
        self.ff = pms.FeatureFinder()
        self.features = pms.FeatureMap()
        self.name = pms.FeatureFinderAlgorithmPicked().getProductName()
        self.params = pms.FeatureFinder().getParameters(self.name)
    
    
    def main(self):
        
        self.peak_picking()
        self.feature_detection()
    
    
    def peak_picking(self):

        pms.MzMLFile().load(self.raw_mzml, self.raw_map)
        pp = pms.PeakPickerHiRes()
        pp.pickExperiment(self.raw_map, self.picked_out_map)
        self.picked_out_map.updateRanges()
        pms.MzMLFile().store(self.peaks_file, self.picked_out_map)
    
    
    def feature_detection(self):
        
        self.options.setMSLevels([1,1])
        self.fh.setOptions(self.options)
        self.fh.load(self.peaks_file, self.picked_input_map)
        self.picked_input_map.updateRanges()
        self.ff.setLogType(pms.LogType.CMD)
        
        # Run the feature finder
        self.ff.run(
            self.name,
            self.picked_input_map,
            self.features,
            self.params,
            self.seeds
        )
        self.features.setUniqueIds()
        self.fh = pms.FeatureXMLFile()
        self.fh.store(self.features_file, self.features)
    
    
    def reload(self):

        modname = self.__class__.__module__
        mod = __import__(modname, fromlist = [modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)


if __name__ == "__main__":
    
    a = Data_Processing(
        raw_mzml = (
            "/home/igor/Documents/Black_scripts/Raw__data/"
            "STARD10_invivo_raw/mzml/"
            "150310_Popeye_MLH_AC_STARD10_A10_pos.mzML"
        ),
        peaks_file = "Test_STARD10_A10_pos_picked_HiRes.mzML",
        features_file = "Test_featureXMLmap.featureXML "
    )
    
    a.peak_picking()
    a.feature_detection()
