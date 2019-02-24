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


class Data_Processing(object):
    """
    Main class for ms data processing using pyopenms library.

    Parameters
    ----------
    
    """
    def __init__(
            self,
            raw_data_file = None,
            peaked_file = None,
            feature_file = None,
            seeds = None,
            fh = None,
            options = None,
            picked_input_map = None,
            ff = None,
            features = None,
            name = None,
            params = None,
            featureXML_file_name=None,
            raw_map = None
        ):
        
        self.featureXML_file_name = featureXML_file_name
        self.raw_data_file = raw_data_file
        self.peaked_file = peaked_file
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
        

    def peak_picking(self):

        pms.MzMLFile().load(self.raw_data_file, self.raw_map)
        pp = pms.PeakPickerHiRes()
        pp.pickExperiment(self.raw_map, self.picked_out_map)
        self.picked_out_map.updateRanges()
        pms.MzMLFile().store(self.peaked_file, self.picked_out_map)

    def feature_detection(self):

        self.options.setMSLevels([1,1])
        self.fh.setOptions(self.options)
        self.fh.load(self.peaked_file, self.picked_input_map)
        self.picked_input_map.updateRanges()
        self.ff.setLogType(pms.LogType.CMD)
        
        # Run the feature finder
        self.ff.run(self.name, self.picked_input_map, self.features, self.params, self. seeds)
        self.features.setUniqueIds()
        self.fh = pms.FeatureXMLFile()
        self.fh.store(self.featureXML_file_name, self.features)
        

    def reload(self):

        modname = self.__class__.__module__
        mod = __import__(modname, fromlist = [modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)


if __name__ == "__main__":
    
    a = Data_Processing(raw_data_file = "/home/igor/Documents/Black_scripts/Raw__data/STARD10_invivo_raw/mzml/150310_Popeye_MLH_AC_STARD10_A10_pos.mzML",
                         peaked_file = "Test_STARD10_A10_pos_picked_HiRes.mzML",
                        featureXML_file_name = "Test_featureXMLmap.featureXML ")
    
    a.peak_picking()
    a.feature_detection()
