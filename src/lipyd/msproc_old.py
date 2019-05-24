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
import lipyd.session as session
import lipyd.common as common
import lipyd.settings as settings


class MSPreprocess(session.Logger):
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
            centroid_input_map = None,
            ff = None,
            features = None,
            name = None,
            params = None,
            profile_map = None,
            smooth_profile_data = False,
            peak_picking_iterative = True,
            feature_min_spectra = 5,
            feature_min_rt_span = .09,
            feature_max_rt_span = 4.5,
            feature_finder_param = None,
            gaussian_smoothing_param = None,
            export_smoothed_profile = False,
            logger = None,
            console = False,
            mass_traces = None,
            mtd_process = None,
            mtd_params = None,
            splitted_mt = None,
            epdet_process = None,
            epdet_params = None,
            chromatograms = None,
            ffm_process = None,
            ffm_params = None,
            feature_map = None,
            reference = None,
            toAlign = None,
            xml_file = None, 
            ma_algorithm = None,
            ma_params = None,
            transformation = None,
            aligned_fm_1 = None,
            aligned_fm_2 = None,
            input_fm_1 = None,
            input_fm_2 = None
        ):
        
        session.Logger.__init__(self, name = 'ms_preprocess')
        
        # filenames and paths
        self.profile_mzml = profile_mzml
        self.centroid_mzml = centroid_mzml
        self.features_file = features_file
        self.wd_root = wd
        self.smooth_profile_data = smooth_profile_data
        self.peak_picking_iterative = peak_picking_iterative
        self.feature_min_spectra = feature_min_spectra
        self.feature_min_rt_span = feature_min_rt_span
        self.feature_max_rt_span = feature_max_rt_span
        self.ff_param = feature_finder_param
        self.gs_param = gaussian_smoothing_param
        self.export_smoothed_profile = export_smoothed_profile
        self.console = console
        self.mass_traces = mass_traces
        self.mtd_process = mtd_process
        self.mtd_params = mtd_params
        self.splitted_mt = splitted_mt,
        self.epdet_process = epdet_process,
        self.epdet_params = epdet_params
        self.chromatograms = chromatograms,
        self.ffm_process = ffm_process,
        self.ffm_params = ffm_params,
        self.feature_map = feature_map

        self.reference = reference
        self.toAlign = toAlign
        self.xml_file = xml_file
        self.ma_algorithm = ma_algorithm
        self.ma_params = ma_params
        self.transformation = transformation
        self.aligned_fm_1 = aligned_fm_1
        self.aligned_fm_2 = aligned_fm_2
        self.input_fm_1 = input_fm_1
        self.input_fm_2 = input_fm_2
    
    
    def main(self):
        
        self.setup()
        self.peak_picking()
        self.feature_finding_metabo()
    
    
    def setup(self):
        
        self.set_paths()
        self.set_feature_finder_param()
    
    
    def set_paths(self):
        
        # default name for all files:
        # name of the input mzML with the path and extension removed
        if not hasattr(self, 'name'):
            
            input_file = self.profile_mzml or self.centroid_mzml
            
            self.name = os.path.splitext(os.path.basename(input_file))[0]
        
        # the working directory
        self.wd_root = self.wd_root or settings.get('ms_preproc_wd')
        self.wd = os.path.join(self.wd_root, self.name)
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
                
                self._log(
                    'Warning: unknown parameter for '
                    'pyopenms.FeatureFinder: `%s`' % param.decode('ascii')
                )
                
                continue
            
            self.ff_param_oms.setValue(param, value)
    
    
    def peak_picking(self):
        
        if not self.profile_mzml:
            
            self._log('No profile data provided, not doing peak picking.')
            
            return
        
        self._log(
            'Performing peak picking on experiment `%s`.' % self.name
        )
        
        self.open_profile_mzml()
        self.smoothing()
        self.run_peak_picker()
        self.export_centroid_mzml()
        
        
        self._log(
            'Peak picking finished. Centroid data has been '
            'written to `%s`.' % self.centroid_mzml
        )
    
    
    def open_profile_mzml(self):
        
        self._log('Loading profile data from `%s`.' % self.profile_mzml)
        
        # As I understand this belongs to the peak picking
        # hence I moved here, we don't need these attributes in __init__
        self.profile_map = oms.MSExperiment()
        oms.MzMLFile().load(self.profile_mzml, self.profile_map)
    
    
    def smoothing(self):
    
        if not self.smooth_profile_data:
            
            return
        
        self._log('Smoothing profile data by Gaussian filter.')
        
        gs = oms.GaussFilter()
        param = gs.getDefaults()
        self._oms_set_param(self.gs_param, param)
        gs.setParameters(param)
        
        gs.filterExperiment(self.profile_map)
        
        if self.export_smoothed_profile:
            
            self.smoothed_profile_mzml = os.path.join(
                self.wd,
                '%s__smoothed.mzML' % self.name,
            )
            self._log(
                'Exporting smoothed profile '
                'data into `%s`.' % self.smoothed_profile_mzml
            )
            oms.MzMLFile().store(self.smoothed_profile_mzml, self.profile_map)
    
    
    def run_peak_picker(self):
        
        self._log('Starting peak picking.')
        
        if self.peak_picking_iterative:
            
            self.pp = oms.PeakPickerIterative()
            
        else:
            
            self.pp = oms.PeakPickerHiRes()
        
        self.centroid_out_map = oms.MSExperiment()
        
        self.pp.pickExperiment(self.profile_map, self.centroid_out_map)
        
        self.centroid_out_map.updateRanges()
    
    
    def export_centroid_mzml(self):
        
        oms.MzMLFile().store(self.centroid_mzml, self.centroid_out_map)
    
    
    def feature_detection(self):
        
        if not self.centroid_mzml:
            
            self._log(
                'No centroid data available, '
                'unable to run feature detection.'
            )
            
            return
        
        self._log('Starting feature detection.')
        
        # As I understand this belongs to the feature detection
        # hence I moved here, we don't need these attributes in __init__
        self.open_centroid_mzml()
        self.setup_feature_finder()
        
        self.run_feature_finder()
        self.export_features()
        
        self._log('Feature detection finished.')
    
    
    def open_centroid_mzml(self):
        
        self._log('Loading centroid data from `%s`.' % self.centroid_mzml)
        
        # opening and reading centroided data from mzML
        self.centroid_mzml_fh = oms.MzMLFile()
        self.centroid_input_map = oms.MSExperiment()
        self.centroid_mzml_options = oms.PeakFileOptions()
        self.centroid_mzml_options.setMaxDataPoolSize(10000)
        self.centroid_mzml_options.setMSLevels([1,1])
        self.centroid_mzml_fh.setOptions(self.centroid_mzml_options)
        self.centroid_mzml_fh.load(
            self.centroid_mzml,
            self.centroid_input_map,
        )
        self.centroid_input_map.updateRanges()
    
    
    def setup_feature_finder(self):
        
        # setting up the FeatureFinder
        self.seeds = oms.FeatureMap()
        self.ff = oms.FeatureFinder()
        self.features = oms.FeatureMap()
        self.ff.setLogType(oms.LogType.CMD)
    
    
    def run_feature_finder(self):
        
        self._log('Running the feature finder.')
        
        # running the feature finder
        self.ff.run(
            self.ff_type,
            self.centroid_input_map,
            self.features,
            self.ff_param_oms,
            self.seeds,
        )
        
        # saving features into featureXML
        self.features.setUniqueIds()
    
    def export_features(self):
        
        self._log('Saving features into `%s`.' % self.features_file)
        self.features_xml = oms.FeatureXMLFile()
        self.features_xml.store(self.features_file, self.feature_map)
    
    def export_chromatograms_data(self):

        self._log('Saving chromatograms into `%s`.' % self.features_file)
        self.chromatogram_mzml = oms.FeatureXMLFile()
        self.chromatogram_mzml.store(self.chromatograms, self.feature_map)


    
    @staticmethod
    def _dict_ensure_bytes(d):
        
        return common.dict_ensure_bytes(d)
    
    
    @staticmethod
    def _oms_set_param(param, target):
        
        param = common.dict_ensure_bytes(param)
        all_param = target.keys()
        
        for par, value in iteritems(param):
            
            target.setValue(par, value)
    

    def map_alignment(self, **kwargs):
        
        self.load_feature_maps()

        self.choose_ma_algorithm()

        self.set_ma_parameters()

        self.run_ma()

        self.store_aligned_maps()
    
    
    def load_feature_maps(self, **kwargs):
        
        self.reference = oms.FeatureMap()
        self.toAlign = oms.FeatureMap()
        self.xml_file = oms.FeatureXMLFile()
        self.xml_file.load(self.input_fm_1, self.reference)
        self.xml_file.load(self.input_fm_2, self.toAlign)
        

    def choose_ma_algorithm(self, **kwargs):

        #create map alignment algorithm
        self.ma_algorithm = oms.MapAlignmentAlgorithmPoseClustering()

    
    def set_ma_parameters(self, **kwargs):
        
        #set parameters
        self.ma_params = self.ma_algorithm.getParameters()    # oms.Param()
        self.ma_params.setValue('superimposer:scaling_bucket_size', 0.01)
        self.ma_params.setValue('superimposer:shift_bucket_size', 0.1)
        #p.setValue('superimposer:num_used_points', 3) #only use first three points -> different results as before expected
        self.ma_params.setValue('superimposer:max_shift', 2000.0 )
        self.ma_params.setValue('superimposer:max_scaling', 2. )
        self.ma_params.setValue('max_num_peaks_considered', -1 ) # -1 is all;
        self.ma_algorithm.setParameters(self.ma_params)
        self.ma_algorithm.setReference(self.reference)
    
    
    def run_ma(self, **kwargs):
        
        #create object for the computed transformation
        self.transformation = oms.TransformationDescription()
        #align
        self.ma_algorithm.align(self.toAlign, self.transformation)
    
    
    def store_aligned_maps(self, **kwargs):
        
        #store results
        #self.xml_file.store(self.aligned_fm_1, self.reference)
        self.xml_file.store(self.aligned_fm_2, self.toAlign)
    
    
    def feature_finding_metabo(self, export_chromatograms=False):
        
        self.load_centroid_mzml()

        self.mt_detection()

        self.elution_prof_detection()

        self.ff_metabo()

        self.export_features()

        if export_chromatograms is True:

            self.export_chromatograms_data()

    def load_centroid_mzml(self):
        
        self.centroid_input_map = oms.PeakMap()
        oms.MzMLFile().load(self.centroid_mzml, self.centroid_input_map)
        self.feature_map = oms.FeatureMap()

    def mt_detection(self):

        self.mass_traces = []

        self.mtd_process = oms.MassTraceDetection()
        
        self.mtd_params = self.mtd_process.getDefaults()
        
        self.mtd_params.setValue('noise_threshold_int', 10.0 )  
        self.mtd_params.setValue('chrom_peak_snr', 3.0 )
        self.mtd_params.setValue('mass_error_ppm', 20.0 )
        self.mtd_params.setValue('reestimate_mt_sd', b'true' )
        self.mtd_params.setValue('quant_method', b'area')
        
        self.mtd_process.setParameters(self.mtd_params)
        
        self.mtd_process.run(self.centroid_input_map, self.mass_traces)
    
    
    def elution_prof_detection(self, **kwargs):
        
        self.splitted_mt = []
        
        self.epdet_process = oms.ElutionPeakDetection()
        
        self.epdet_params = self.epdet_process.getDefaults()
        
        self.epdet_params.setValue('chrom_peak_snr', 3.0 )
        self.epdet_params.setValue('chrom_fwhm', 5.0 )
        self.epdet_params.setValue('width_filetering', b'fixed')
        
        self.epdet_process.setParameters(self.epdet_params)
        
        self.epdet_process.detectPeaks(self.mass_traces, self.splitted_mt)
        
        
    def ff_metabo(self):
        
        self.chromatograms = [[]]
        
        self.ffm_process = oms.FeatureFindingMetabo()
        
        self.ffm_params = self.ffm_process.getDefaults()
        
        self.ffm_params.setValue('charge_lower_bound', 1 )
        self.ffm_params.setValue('charge_upper_bound', 1 )        
        self.ffm_params.setValue('enable_RT_filtering', b'true') 
        self.ffm_params.setValue('isotope_filtering_model', b'metabolites (5% RMS)') 
        self.ffm_params.setValue('mz_scoring_13C', b'true') 
        self.ffm_params.setValue('report_convex_hulls', b'true') 
        self.ffm_params.setValue('remove_single_traces', b'true') 
        
        self.ffm_process.setParameters(self.ffm_params)
        
        self.ffm_process.run(self.splitted_mt, self.feature_map, self.chromatograms)
    
    
    def reload(self):

        modname = self.__class__.__module__
        mod = __import__(modname, fromlist = [modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)


if __name__ == "__main__":
    
    a = MSPreprocess(
        profile_mzml = (
            "/home/igor/Documents/Scripts/Raw_data_STARD10/"
           "150310_Popeye_MLH_AC_STARD10_A10_pos.mzML"),
        centroid_mzml = "/home/igor/Documents/lipyd/src/lipyd_ms_preproc/150310_Popeye_MLH_AC_STARD10_A09_pos/150310_Popeye_MLH_AC_STARD10_A10_pos__peaks.mzML",
        input_fm_1 = "/home/igor/Documents/lipyd/src/lipyd_ms_preproc/150310_Popeye_MLH_AC_STARD10_A09_pos/150310_Popeye_MLH_AC_STARD10_A09_pos__features.featureXML", 
        input_fm_2 = "/home/igor/Documents/lipyd/src/lipyd_ms_preproc/150310_Popeye_MLH_AC_STARD10_A09_pos/150310_Popeye_MLH_AC_STARD10_A10_pos__features.featureXML",
        aligned_fm_1 = "/home/igor/Documents/lipyd/src/lipyd_ms_preproc/150310_Popeye_MLH_AC_STARD10_A09_pos/test_aligned_map_A09.featureXML",
        aligned_fm_2 = "/home/igor/Documents/lipyd/src/lipyd_ms_preproc/150310_Popeye_MLH_AC_STARD10_A10_pos/test_aligned_map_A10.featureXML" 
    )
    a.setup()
    a.peak_picking()
    a.feature_finding_metabo()
    a.map_alignment()
