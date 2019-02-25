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

import os
import imp

import pyopenms as oms

#TODO: Please Igor integrate logging
# At each major step send messages to the logger
# As an example see
# https://git.embl.de/grp-gavin/ltp_ms/blob/ltp/ltp/featureproc.py
import lipyd.log as log
import lipyd.settings as settings


class MSPreprocess(object):
    #TODO: Igor, please describe the parameters, I don't understand half
    # of them :)
    # Also for many I don't see where we use them.
    
    """
    Main class for ms data processing using pyopenms library.

    Parameters
    ----------
    """
    
    
    def __init__(
            self,
            profile_mzml,
            peaks_file = None,
            features_file = None,
            wd = None,
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
        
        # filenames and paths
        self.profile_mzml = profile_mzml
        self.peaks_file = peaks_file
        self.features_file = features_file
        self.wd = wd
        
        self.set_paths()
    
    
    def set_paths(self):
        
        # default name for all files:
        # name of the input mzML with the path and extension removed
        self._name = '.'.join(
            os.path.split(self.profile_mzml)[-1].split('.')[:-1]
        )
        
        # the working directory
        self.wd = self.wd or settings.get('ms_preproc_wd')
        self.wd = os.path.join(self.wd, self._name)
        os.makedirs(self.wd, exist_ok = True)
        
        self.peaks_file = self.peaks_file or '%s__peaks.mzML' % self._name
        self.peaks_file = os.path.join(wd)
        
        self.features_file = (
            self.features_file or '%s__features.featureXML' % self._name
        )
    
    
    def main(self):
        
        self.peak_picking()
        self.feature_detection()
    
    
    def peak_picking(self):
        
        # As I understand this belongs to the peak picking
        # hence I moved here, we don't need these attributes in __init__
        self.raw_map = oms.MSExperiment()
        self.picked_out_map = oms.MSExperiment()
        
        oms.MzMLFile().load(self.profile_mzml, self.raw_map)
        pp = oms.PeakPickerHiRes()
        pp.pickExperiment(self.raw_map, self.picked_out_map)
        self.picked_out_map.updateRanges()
        oms.MzMLFile().store(self.peaks_file, self.picked_out_map)
    
    
    def feature_detection(self):
        
        # As I understand this belongs to the feature detection
        # hence I moved here, we don't need these attributes in __init__
        self.seeds = oms.FeatureMap()
        self.fh = oms.MzMLFile()
        self.options = oms.PeakFileOptions()
        self.picked_input_map = oms.MSExperiment()
        self.ff = oms.FeatureFinder()
        self.features = oms.FeatureMap()
        self.name = oms.FeatureFinderAlgorithmPicked().getProductName()
        self.params = oms.FeatureFinder().getParameters(self.name)
        
        self.options.setMSLevels([1,1])
        self.fh.setOptions(self.options)
        self.fh.load(self.peaks_file, self.picked_input_map)
        self.picked_input_map.updateRanges()
        self.ff.setLogType(oms.LogType.CMD)
        
        # Run the feature finder
        self.ff.run(
            self.name,
            self.picked_input_map,
            self.features,
            self.params,
            self.seeds,
        )
        self.features.setUniqueIds()
        self.fh = oms.FeatureXMLFile()
        self.fh.store(self.features_file, self.features)
    
    
    def reload(self):

        modname = self.__class__.__module__
        mod = __import__(modname, fromlist = [modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)


if __name__ == "__main__":
    
    a = Data_Processing(
        profile_mzml = (
            "/home/igor/Documents/Black_scripts/Raw__data/"
            "STARD10_invivo_raw/mzml/"
            "150310_Popeye_MLH_AC_STARD10_A10_pos.mzML"
        ),
        peaks_file = "Test_STARD10_A10_pos_picked_HiRes.mzML",
        features_file = "Test_featureXMLmap.featureXML "
    )
    
    a.peak_picking()
    a.feature_detection()
