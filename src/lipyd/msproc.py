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

from future.utils import iteritems

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
            profile_mzml = None,
            centroid_mzml = None,
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
            raw_map = None,
            smooth_profile_data = False,
            peak_picking_iterative = True,
            feature_min_spectra = 5,
            feature_min_rt_span = .9,
            feature_max_rt_span = 4.5,
            feature_finder_param = None,
            gaussian_smoothing_param = None,
            logger = None,
        ):
        
        self.log = (
            logger or
            log.new_logger(
                'lipyd.msproc',
                logdir = 'lipyd_log',
            )
        )
        
        # filenames and paths
        self.profile_mzml = profile_mzml
        self.centroid_mzml = centroid_mzml
        self.features_file = features_file
        self.wd = wd
        self.smooth_profile_data = smooth_profile_data
        self.peak_picking_iterative = peak_picking_iterative
        self.feature_min_spectra = feature_min_spectra
        self.feature_min_rt_span = feature_min_rt_span
        self.feature_max_rt_span = feature_max_rt_span
        self.ff_param = feature_finder_param
        self.gs_param = gaussian_smoothing_param
    
    
    def main(self):
        
        self.setup()
        self.peak_picking()
        self.feature_detection()
    
    
    def setup(self):
        
        self.set_paths()
        self.set_feature_finder_param()
    
    
    def set_paths(self):
        
        # default name for all files:
        # name of the input mzML with the path and extension removed
        if not hasattr(self, 'name'):
            
            input_file = self.profile_mzml or self.centroid_mzml
            
            self.name = '.'.join(
                os.path.splitext(os.path.basename(input_file))[0]
            )
        
        # the working directory
        self.wd = self.wd or settings.get('ms_preproc_wd')
        self.wd = os.path.join(self.wd, self.name)
        os.makedirs(self.wd, exist_ok = True)
        
        self.centroid_mzml = self.centroid_mzml or '%s__peaks.mzML' % self.name
        self.centroid_mzml = os.path.join(self.wd, self.centroid_mzml)
        
        self.features_file = (
            self.features_file or '%s__features.featureXML' % self.name
        )
        self.features_file = os.path.join(self.wd, self.features_file)
    
    
    def set_feature_finder_param(self):
        
        self.ff_param = self.ff_param or {}
        
        self.ff_param = dict(
            (
                param.encode('ascii') if hasattr(param, 'encode') else param,
                value,
            )
            for param, value in iteritems(self.ff_param)
        )
        
        extra_param = dict((
            (b'mass_trace:min_spectra', self.feature_min_spectra),
            (b'feature:min_rt_span', self.feature_min_rt_span),
            (b'feature:max_rt_span', self.feature_max_rt_span),
        ))
        
        extra_param.update(self.ff_param)
        
        self.ff_param = extra_param
        
        self.ff_type = oms.FeatureFinderAlgorithmPicked().getProductName()
        self.ff_param_oms = oms.FeatureFinder().getParameters(self.ff_type)
        all_param = self.ff_param_oms.keys()
        
        for param, value in iteritems(self.ff_param):
            
            if param not in all_param:
                
                self.log.msg(
                    'Warning: unknown parameter for '
                    'pyopenms.FeatureFinder: `%s`' % param.decode('ascii')
                )
                
                continue
            
            self.ff_param_oms.setValue(param, value)
    
    
    def peak_picking(self):
        
        if not self.profile_mzml:
            
            self.log.msg('No profile data provided, not doing peak picking.')
            
            return
        
        self.log.msg(
            'Performing peak picking on experiment `%s`.' % self.name
        )
        self.log.msg('Using profile data from `%s`.' % self.profile_mzml)
        self.log.msg(
            'Centroid data will be written to `%s`.' % self.centroid_mzml
        )
        
        # As I understand this belongs to the peak picking
        # hence I moved here, we don't need these attributes in __init__
        self.raw_map = oms.MSExperiment()
        self.picked_out_map = oms.MSExperiment()
        
        oms.MzMLFile().load(self.profile_mzml, self.raw_map)
        
        if self.smooth_profile_data:
            
            self.log.msg('Smoothing profile data by Gaussian filter.')
            
            gs = oms.GaussFilter()
            gs.filterExperiment(self.raw_map)
        
        self.log.msg('Starging peak picking.')
        
        if self.peak_picking_iterative:
            
            self.pp = oms.PeakPickerIterative()
            
        else:
            
            self.pp = oms.PeakPickerHiRes()
        
        self.pp.pickExperiment(self.raw_map, self.picked_out_map)
        self.picked_out_map.updateRanges()
        oms.MzMLFile().store(self.centroid_mzml, self.picked_out_map)
        
        self.log.msg(
            'Peak picking finished. Centroid data has been '
            'written to `%s`.' % self.centroid_mzml
        )
    
    
    def feature_detection(self):
        
        self.log.msg('Starting feature detection.')
        
        # As I understand this belongs to the feature detection
        # hence I moved here, we don't need these attributes in __init__
        self.open_centroid_mzml()
        self.setup_feature_finder()
        
        self.run_feature_finder()
        self.save_features()
        
        self.log.msg('Feature detection finished.')
    
    
    def open_centroid_mzml(self):
        
        self.log.msg('Loading centroid data from `%s`.' % self.centroid_mzml)
        
        # opening and reading centroided data from mzML
        self.centroid_mzml_fh = oms.MzMLFile()
        self.centroid_mzml_options = oms.PeakFileOptions()
        self.centroid_mzml_options.setMaxDataPoolSize(10000)
        self.centroid_mzml_options.setMSLevels([1,1])
        self.fh.setOptions(self.centroid_mzml_options)
        self.centroid_mzml_fh.load(self.centroid_mzml, self.picked_input_map)
        self.picked_input_map.updateRanges()
    
    
    def setup_feature_finder(self):
        
        # setting up the FeatureFinder
        self.seeds = oms.FeatureMap()
        self.picked_input_map = oms.MSExperiment()
        self.ff = oms.FeatureFinder()
        self.features = oms.FeatureMap()
        self.ff.setLogType(oms.LogType.CMD)
    
    
    def run_feature_finder(self):
        
        self.log.msg('Running the feature finder.')
        
        # running the feature finder
        self.ff.run(
            self.ff_type,
            self.picked_input_map,
            self.features,
            self.ff_param_oms,
            self.seeds,
        )
        
        # saving features into featureXML
        self.features.setUniqueIds()
    
    def save_features(self):
        
        self.log.msg('Saving features into `%s`.' % self.features_file)
        
        # Igor, here you had self.fh, file handler for the featureXML
        # 10 lines above same variable name was the mzML file
        # very confusing to use the same name for different things! :)
        self.features_xml_fh = oms.FeatureXMLFile()
        self.features_xml_fh.store(self.features_file, self.features)
    
    
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
        centroid_mzml = "Test_STARD10_A10_pos_picked_HiRes.mzML",
        features_file = "Test_featureXMLmap.featureXML "
    )
    
    a.peak_picking()
    a.feature_detection()
