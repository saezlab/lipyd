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
#  Website: https://saezlab.github.io/lipyd
#

from future.utils import iteritems

import os
import sys
import imp
import re
import copy

import numpy as np
import pyopenms as oms

import lipyd.session as session
import lipyd.common as common
import lipyd.settings as settings

OPENMS_OBJ_TYPES = (
    oms.PeakMap,
    oms.FeatureMap,
    oms.MSExperiment,
    oms.ConsensusMap,
)


class MethodParamHandler(session.Logger):
    """
    Base class for setup and check parameters for all derived classes.
    
    Parameters
    ----------
    **kwargs
        Parameters for the OpenMS object as keyword arguments.

    Attributes
    ----------
    openms_method : class
        OpenMS class having ``setParameters`` method.
    openms_obj : object
        OpenMS object having ``setParameters`` method.
    """
    
    
    _module_param_keys = {
        oms.PeakPickerHiRes: (
            'peak_picking_param',
        ),
        oms.ElutionPeakDetection: (
            'feature_finding_common_param',
            'elution_peak_detection_param',
        ),
        oms.MassTraceDetection: (
            'feature_finding_common_param',
            'mass_trace_detection_param',
        ),
        oms.FeatureFindingMetabo: (
            'feature_finding_common_param',
            'feature_finding_metabo_param',
        ),
        oms.MapAlignmentAlgorithmPoseClustering: (
            'map_alignment_param',
        ),
    }
    
    
    def __init__(
            self,
            openms_method,
            name = 'param_handler',
            **kwargs
        ):
        
        # most of the times log initialized in the OpenmsMethodWrapper
        if not hasattr(self, '_log_name'):
            
            session.Logger.__init__(self, name = name)
        
        self.openms_method = openms_method
        self._kwargs_param = kwargs
    
    
    def main(self):
        
        self.setup()
    
    
    def setup(self):
        """
        Collects and updates parameters according to the OpenMS object type.
        """
        
        self.get_openms_object_name()
        self.collect_param()
        self.set_param()
    
    
    def get_openms_object_name(self):
        """
        Creates a string representation of the OpenMS object type.
        """
        
        self.openms_object_type = '%s.%s' % (
            self.openms_method.__module__,
            self.openms_method.__name__,
        )
    
    
    def collect_param(self):
        """
        Collects parameters from (1) keyword arguments for this class,
        (2) settings of the `lipyd` module and (3) OpenMS defaults in this
        order of priority.
        """
        
        self.get_openms_defaults()
        self.get_module_param()
        
        _kwargs_param = self.process_param_dict(self._kwargs_param)
        _module_param = self.process_param_dict(self._module_param)
        self.param = {}
        
        self.param.update(self._openms_defaults)
        self.param.update(_module_param)
        self.param.update(_kwargs_param)
    
    
    def set_param(self):
        """
        Sets the parameters of the OpenMS object according to the current
        `param` dictionary.
        """
        
        self.openms_param = self.openms_obj.getDefaults()
        
        for key, value in iteritems(self.param):
            
            self.openms_param.setValue(key, value)
            
            self._log(
                'Parameter `{}` of `{}` set to `{}`.'.format(
                    common.ensure_unicode(key),
                    self.openms_object_type,
                    common.ensure_unicode(value),
                ),
                1,
            )
        
        self.openms_obj.setParameters(self.openms_param)
    
    
    @staticmethod
    def bytes_or_numeric(value, float_ok = True):
        """
        Returns the `value` in a bytes string or a numeric representation.
        """
        
        if isinstance(value, bool):
            
            return value.__str__().lower().encode()
            
        elif isinstance(value, np.float64):
            
            return float(value)
            
        elif isinstance(value, (int, float)):
            
            return value
            
        elif (common.is_numeric(value) and float_ok) or common.is_int(value):
            
            return common.to_number(value)
            
        elif hasattr(value, 'encode'):
            
            return value.encode()
            
        else:
            
            return value
    
    
    def process_param_dict(self, param):
        """
        Makes sure all keys and values in the dict are byte strings or
        numeric.
        """
        
        return dict(
            (
                self.bytes_or_numeric(k, float_ok = False),
                self.bytes_or_numeric(v)
            )
            for k, v in iteritems(param)
        )
    
    
    def get_openms_defaults(self):
        """
        Queries the default parameter values from the OpenMS object and
        converts them to a `dict`.
        """
        
        # get the defaults from openms
        self._openms_defaults = dict(self.openms_obj.getDefaults())
    
    
    def get_module_param(self):
        """
        Queries the current module level settings of `lipyd` for the
        particular OpenMS class.
        """
        
        self._module_param = {}
        
        if self.openms_method in self._module_param_keys:
            
            for key in self._module_param_keys[self.openms_method]:
                
                self._module_param.update(settings.get(key))
            
        else:
            
            self._log(
                'Could not find settings key for OpenMS '
                'object type `%s`.' % self.openms_object_type
            )


class PathHandlerBase(session.Logger):
    
    _method_keys = {
        oms.PeakPickerHiRes: 'centroided',
        oms.FeatureFindingMetabo: 'feature',
        oms.MapAlignmentAlgorithmPoseClustering: 'aligned',
    }
    
    def __init__(
            self,
            input_path = None,
            input_obj = None,
            method_key = None,
            sample_id = None,
            sample_id_method = None,
            input_ext = None,
            input_filter = None,
            name = 'path_handler',
            n_experiment = None,
            n_sample = None,
            multi_file_input = False,
        ):
        
        # most of the times log initialized in the OpenmsMethodWrapper
        if not hasattr(self, '_log_name'):
            
            session.Logger.__init__(self, name = name)
        
        self.multi_file_input = multi_file_input
        self._input_path = input_path
        self.input_obj = input_obj
        self.method_key = method_key
        self._input_ext = input_ext
        self._input_filter = input_filter
        self._sample_id_method = sample_id_method
        self.sample_id = sample_id
        self.n_experiment = n_experiment
        self.n_sample = n_sample
    
    
    def main(self):
        
        self.set_paths()
    
    
    def set_paths(self):
        
        if self._input_path is not None:
            
            self._set_input_path()
        
        self._set_input()
    
    
    def _set_input_path(self):
        
        if isinstance(self._input_path, common.basestring):
            
            # check if this is a directory
            if os.path.isdir(self._input_path):
                
                # look up all files in the directory
                self._input_path = [
                    os.path.join(self._input_path, fname)
                    for fname in os.listdir(self._input_path)
                ]
        
        # if we have something else than a single existing file path
        self.multi_file_input = (
            isinstance(self._input_path, (tuple, list, np.ndarray))
        )
        
        if self.multi_file_input:
            
            self._filter_input_paths()
            # the `input_path` will be just the first file
            # just to make it easier for the downstream methods
            self.input_path = self._input_path[0]
            
        else:
            
            # the simple case when we have only one file
            self.input_path = self._input_path
    
    
    def _filter_input_paths(self):
        
        if isinstance(self._input_filter, common.basestring):
            
            refname = re.compile(self._input_filter)
            
            def _input_filter(path):
                
                return bool(refname.fullmatch(path))
            
            self._input_filter = _input_filter
        
        self._input_path = [
            path
            for path in self._input_path
            if (
                self._input_ext is None or
                os.path.splitext(path)[1][1:] == self._input_ext
            ) and (
                self._input_filter is None or
                # this method recieves only the file name
                self._input_filter(os.path.basename(path))
            )
        ]
    
    
    def _set_input(self):
        
        if (
            isinstance(self.input_obj, OPENMS_OBJ_TYPES) or
            (
                isinstance(self.input_obj, (tuple, list, np.ndarray)) and
                all(
                    isinstance(obj, OPENMS_OBJ_TYPES)
                    for obj in self.input_obj
                )
            )
        ):
            
            self._input = self.input_obj
            
        elif (
            isinstance(self.input_path, common.basestring) or
            (
                isinstance(self.input_path, (tuple, list, np.ndarray)) and
                all(
                    isinstance(obj, common.basestring)
                    for path in self.input_path
                )
            )
        ):
            
            self._input = self.input_path
    
    
    def _set_sample_id(self):
        
        if not self.sample_id:
            
            if self.multi_file_input:
                
                self.sample_id = [
                    self._get_sample_id(path)
                    for path in self._input_path
                ]
                
            else:
                
                self.sample_id = self._get_sample_id(self._input_path)
    
    
    def _get_sample_id(self, _input):
        """
        Parameters
        ----------
        _input : str,object
            Either the path to the input file or an OpenMS object.
        """
        
        name = None
        sample_id = None
        
        # first we try to get a string which carries information
        # regarding the identity of the sample
        if isinstance(self._sample_id_method, common.basestring):
            
            sample_id = self._sample_id_method
            
        else:
            
            if isinstance(_input, OPENMS_OBJ_TYPES):
                
                if hasattr(_input, 'getSourceFiles'):
                    
                    files = _input.getSourceFiles()
                    
                    if files:
                        
                        name = files[0].getNameOfFile().decode()
                
                if not name and hasattr(_input, 'getSample'):
                    
                    sample = _input.getSample()
                    name = '%s__%s' % (
                        sample.getName().decode(),
                        sample.getNumber().decode(),
                    )
                
            elif (
                isinstance(_input, common.basestring) and
                os.path.exists(_input)
            ):
                
                name = os.path.splitext(os.path.basename(_input))[0]
            
            # once we have a string call the sample_id_method on it
            if isinstance(name, common.basestring):
                
                sample_id = (
                    self._sample_id_method(name)
                        if callable(self._sample_id_method) else
                    name
                )
                
            elif isinstance(self._sample_id_method, common.basestring):
                
                sample_id = self._sample_id_method
        
        if not sample_id and hasattr(self, '_set_output_dir'):
            
            experiment = 1
            sample = 1
            if not hasattr(self, 'method_key'):
                self._set_method_key()
            suffix = settings.get('%s_suffix' % self.method_key) or ''
            rename = re.compile(
                r'experiment(\d+)_sample(\d+)' + suffix + r'\.\w+'
            )
            if not hasattr(self, 'output_dir'):
                self._set_output_dir()
            all_files = os.listdir(self.output_dir)
            files = sorted([f for f in all_files if rename.match(f)])
            
            if files:
                
                last_file = files[-1]
                m = rename.match(last_file).groups()
                experiment = int(m[0])
                sample = int(m[1])
            
            sample_id = 'experiment%03u_sample%03u' % (experiment, sample)
        
        return sample_id


class MethodPathHandler(PathHandlerBase):
    """
    Parameters
    ----------
    input_path : str
        The path of the input file must be provided, otherwise this class
        will do nothing. In case of multiple input files you can provide
        a list of paths or the directory containing the files.
    method_key : str,class
        Either an OpenMS class or a string label. Used to identify the
        method and look up the corresponding settings.
    sample_id_method : callable
        A method which generates sample identifier(s) from the input path(s).
        If not provided, the file name without the extension will be used.
    input_ext : str
        In case of multiple input files only the files with this extension
        (type) will be used.
    input_filter : callable,str
        A method which decides if an input file should be used. Should return
        ``True`` if the file is to be used. Alternatively a regex as a string
        which matches the desired input files.
    
    Attributes
    ----------
    multi_file_input : bool
        Tells if the object has single or multiple input paths.
    """
    
    _outexts = {
        'centroided': 'mzML',
        'feature': 'featureXML',
        'aligned': 'featureXML',
        'msproc': 'featureXML',
        'mgf_export': 'mgf',
    }
    _subdirs = {
        settings.get('%s_dir' % key)
        for key in ('centroided', 'feature', 'aligned', 'mgf_export')
    }
    
    def __init__(
            self,
            input_path = None,
            input_obj = None,
            method_key = None,
            output_path = None,
            output_dir = None,
            sample_id = None,
            sample_id_method = None,
            input_ext = None,
            input_filter = None,
            name = 'path_handler',
            n_experiment = None,
            n_sample = None,
        ):
        
        # most of the times log initialized in the OpenmsMethodWrapper
        if not hasattr(self, '_log_name'):
            
            session.Logger.__init__(self, name = name)
        
        PathHandlerBase.__init__(
            self,
            input_path = input_path,
            input_obj = input_obj,
            method_key = method_key,
            sample_id = sample_id,
            sample_id_method = sample_id_method,
            input_ext = input_ext,
            input_filter = input_filter,
            name = name,
            n_experiment = n_experiment,
            n_sample = n_sample,
        )
        
        self.output_path = output_path
        self.method_key = method_key
        self.output_dir = output_dir
    
    
    def main(self):
        
        self.set_paths()
    
    
    def set_paths(self):
        
        PathHandlerBase.set_paths(self)
        self._set_sample_id()
        
        if not self.output_path:
            
            self._set_method_key()
            self._set_output_dir()
            self._set_output_path()
        
        self._create_output_dir()
        self._tell_paths()
        self._check_paths()
    
    
    def _set_output_dir(self):
        
        output_root = (
            # 1: output path explicitely set
            settings.get('output_path_root') or (
                # 2: it's set to be in the input dir
                os.path.dirname(self.input_path)
                    if (
                        settings.get('output_to_input_dir') or
                        not hasattr(self, 'input_path')
                    ) else
                # 3: by default in current wd
                os.getcwd('output_to_input_dir')
            )
        )
        
        path = os.path.split(output_root)
        
        lipyd_wd = settings.get('lipyd_wd') or ''
        
        if path and path[-1] in self._subdirs:
            
            path = path[:-1]
            path = os.path.split(path[0])
        
        if path and path[-1] == lipyd_wd:
            
            path = path[:-1]
            path = os.path.split(path[0])
        
        method_wd = settings.get('%s_dir' % self.method_key) or ''
        
        path = path + (lipyd_wd, method_wd)
        
        self.output_dir = os.path.join(*path)
    
    
    def _create_output_dir(self):
        
        os.makedirs(self.output_dir, exist_ok = True)
    
    
    def _set_output_path(self):
        
        suffix = settings.get('%s_suffix' % self.method_key) or ''
        
        if isinstance(self.input_path, common.basestring):
            
            self._input_basename, ext = (
                os.path.splitext(os.path.basename(self.input_path))
            )
            
        else:
            
            self._input_basename = self.sample_id
        
        outext = (
            self._outexts[self.method_key]
                if self.method_key in self._outexts else
            'mzML'
        )
        self.outfile = '%s%s.%s' % (self._input_basename, suffix, outext)
        
        self.output_path = os.path.join(self.output_dir, self.outfile)
    
    
    def _set_method_key(self):
        
        if (
            not self.method_key and
            hasattr(self, 'openms_method') and
            self.openms_method in self._method_keys
        ):
            
            self.method_key = self._method_keys[self.openms_method]
    
    
    def _tell_paths(self):
        """
        Writes the paths to the log.
        """
        
        self._log('Reading input from `%s`.' % self.input_path)
        self._log('Will write output to `%s`.' % self.output_path)
    
    
    def _check_paths(self):
        """
        Raising exception if input and output paths are identical.
        """
        
        if self.input_path == self.output_path:
            
            self._log(
                'The input and output files are identical: `%s`. '
                'You can set different output path either directly '
                'providing the `output_path` argument or changing '
                'parameters by the `lipyd.settings` module.' % (
                    self.output_path
                )
            )
            
            raise RuntimeError('Identical in/out path, please check the log.')


class OpenmsMethodWrapper(MethodParamHandler, MethodPathHandler):
    """
    Generic base class for OpenMS method wrappers.
    
    Parameters
    ----------
    method : class
        The OpenMS class.
    infile : str
        The input file must be provided, otherwise no input and output
        assumed and only the parameters will be set.
    name : str
        A name for this method, will appear in the log.
    outfile : str
        Output path provided directly. Otherwise will be set according to
        the current settings.
    method_key : str
        
    """
    
    def __init__(
            self,
            method,
            input_path = None,
            input_obj = None,
            name = None,
            output_path = None,
            method_key = None,
            sample_id = None,
            sample_id_method = None,
            **kwargs
        ):
        
        session.Logger.__init__(self, name = name or 'openms_wrapper')
        
        self.openms_method = method
        
        MethodParamHandler.__init__(
            self,
            openms_method = method,
            **kwargs
        )
        
        MethodPathHandler.__init__(
            self,
            input_path = input_path,
            input_obj = input_obj,
            method_key = method_key,
            output_path = output_path,
            sample_id = sample_id,
            sample_id_method = sample_id_method,
        )
    
    
    def main(self):
        
        self.openms_wrapper_setup()
        self.run()
    
    
    def openms_wrapper_setup(self):
        
        self.create_instance()
        self.setup()
        self.set_paths()
    
    
    def create_instance(self):
        
        self.openms_obj = self.openms_method()
        self._log(
            '`%s.%s` instance created.' % (
                self.openms_method.__module__,
                self.openms_method.__name__,
            )
        )
    
    
    def run(self):
        
        pass


class PeakPickerHiRes(OpenmsMethodWrapper):
    """
    Wrapper around ``pyopenms.PeakPickerHiRes``.
    
    This class implements a fast peak-picking algorithm best suited for high
    resolution MS data (FT-ICR-MS, Orbitrap). In high resolution data, the
    signals of ions with similar mass-to-charge ratios (m/z) exhibit little
    or no overlapping and therefore allow for a clear separation. Furthermore,
    ion signals tend to show well-defined peak shapes with narrow peak width.

    This peak-picking algorithm detects ion signals in profile data and
    reconstructs the corresponding peak shape by cubic spline interpolation.
    Signal detection depends on the signal-to-noise ratio which is adjustable
    by the user (see parameter signal_to_noise). A picked peak's m/z and
    intensity value is given by the maximum of the underlying peak spline.

    So far, this peak picker was mainly tested on high resolution data. With
    appropriate preprocessing steps (e.g. noise reduction and baseline
    subtraction), it might be also applied to low resolution data.

    Parameters
    ----------
    infile : str
        Input file mzML with profile data.
    outfile : str
        The output file can be set directly this way. Alternatively see
        the built in path handling parameters in `settings`.
    **kwargs
        Settings passed directly to ``pyopenms.PeakPickerHiRes``.
    """
    
    
    def __init__(
            self,
            input_path = None,
            input_obj = None,
            output_path = None,
            sample_id = None,
            sample_id_method = None,
            **kwargs,
        ):
        
        OpenmsMethodWrapper.__init__(
            self,
            method = oms.PeakPickerHiRes,
            input_path = input_path,
            input_obj = input_obj,
            output_path = output_path,
            name = 'peak_picker',
            sample_id = sample_id,
            sample_id_method = sample_id_method,
            **kwargs
        )
    
    
    def run(self):
        
        self.read()
        self.pick()
        self.write()
    
    
    def read(self):
        
        if isinstance(self.input_obj, oms.MSExperiment):
            
            self.input_map = self.input_obj
            self._log('Input map provided as `pyopenms.MSExperiment` object.')
            
        else:
            
            self.input_map = oms.MSExperiment()
            oms.MzMLFile().load(self.input_path, self.input_map)
            self._log('Input map created from `%s`.' % self.input_path)
    
    
    def pick(self):
        
        self._log('Starting peak picking.')
        self.output_map = oms.MSExperiment()
        self.openms_obj.pickExperiment(self.input_map, self.output_map)
        self._log('Peak picking ready.')
    
    
    def write(self):
        
        self.output_map.updateRanges()
        oms.MzMLFile().store(self.output_path, self.output_map)
        self._log(
            'Centroided data has been written to `%s`.' % self.output_path
        )


class MassTraceDetection(OpenmsMethodWrapper):
    """
    Wrapper around ``pyopenms.MassTraceDetection``.
    """
    
    
    def __init__(
            self,
            input_map,
            mass_traces = None,
            **kwargs
        ):
        
        self.input_map = input_map
        self.mass_traces = mass_traces if mass_traces is not None else []
        
        OpenmsMethodWrapper.__init__(
            self,
            method = oms.MassTraceDetection,
            name = 'mass_trace_detection',
            **kwargs
        )
    
    
    def main(self):
        
        self.create_instance()
        self.run()
    
    
    def run(self):
        
        self.openms_obj.run(self.input_map, self.mass_traces)


class ElutionPeakDetection(OpenmsMethodWrapper):
    """
    Wrapper around ``pyopenms.ElutionPeakDetection``.
    """
    
    
    def __init__(
            self,
            mass_traces,
            mass_traces_split = None,
            **kwargs
        ):
        
        self.mass_traces = mass_traces
        self.mass_traces_split = (
            mass_traces_split if mass_traces_split is not None else []
        )
        
        OpenmsMethodWrapper.__init__(
            self,
            method = oms.ElutionPeakDetection,
            name = 'elution_peak_detection',
            **kwargs
        )
    
    
    def main(self):
        
        self.create_instance()
        self.run()
    
    
    def run(self):
        
        self.openms_obj.detectPeaks(self.mass_traces, self.mass_traces_split)


class FeatureFindingMetabo(OpenmsMethodWrapper):
    """
    Wrapper around ``pyopenms.FeatureFindingMetabo()``.
    Assembles mass traces showing similar isotope- and elution patterns.
    
    Parameters
    ----------
    input_file : str
        Path to an mzML file with centroided data.
    input_map : pyopenms.PeakMap
        A``pyopenms.PeakMap`` object. If provided, the data will be used
        from this instead of reading from ``input_file``. In this case
        ``input_file`` can be even ``None`` if ``output_file`` provided.
        Otherwise ``input_file`` will be used to create the ``output_file``
        path and name.
    output_file : str
        The output file where the features will be saved in ``featureXML``
        format. If not provided it will be set according to the current
        settings.
    mass_trace_detection_param : dict
        Parameters directly for the ``pyopenms.MassTraceDetection`` class.
    elution_peak_detection_param : dict
        Parameters directly for the ``pyopenms.ElutionPeakDetection`` class.
    feature_finding_metabo_param : dict
        Parameters directly for the ``pyopenms.FeatureFindingMetabo`` class.
    **kwargs
        Common parameters will be passed to all of the OpenMS methods above.
    """
    
    
    def __init__(
        self,
        input_path = None,
        input_obj = None,
        output_path = None,
        sample_id = None,
        sample_id_method = None,
        mass_trace_detection_param = None,
        elution_peak_detection_param = None,
        feature_finding_metabo_param = None,
        **kwargs
    ):
        
        self.mass_trace_detection_param = (
            self._combine_param(kwargs, mass_trace_detection_param)
        )
        self.elution_peak_detection_param = (
            self._combine_param(kwargs, elution_peak_detection_param)
        )
        self.feature_finding_metabo_param = (
            self._combine_param(kwargs, feature_finding_metabo_param)
        )
        
        OpenmsMethodWrapper.__init__(
            self,
            method = oms.FeatureFindingMetabo,
            input_path = input_path,
            input_obj = input_obj,
            output_path = output_path,
            name = 'feature_finding_metabo',
            sample_id = sample_id,
            sample_id_method = sample_id_method,
            **self.feature_finding_metabo_param,
        )
    
    
    def run(self):
        
        self.read()
        self.find_features()
        self.write()
    
    
    def read(self):
        """
        Reads the centroided data.
        """
        
        # first try if OpenMS object provided, a PeakMap:
        if isinstance(self.input_obj, oms.PeakMap):
            
            self.input_map = self.input_obj
            self._log(
                'Centroided data provided as `pyopenms.PeakMap` object.'
            )
            
        else:
            
            self.input_map = oms.PeakMap()
            
            # or an MSExperiment object:
            if isinstance(self.input_obj, oms.MSExperiment):
                
                for chromatogram in self.input_obj.getChromatograms():
                    
                    self.input_map.addChromatogram(chromatogram)
                
                for spectrum in self.input_obj.getSpectra():
                    
                    self.input_map.addSpectrum(spectrum)
                
                self.input_map.updateRanges()
                
                self._log(
                    'Centroided data provided as '
                    '`pyopenms.MSExperiment` object, converted to '
                    '`pyopenms.PeakMap` object.'
                )
                
            # here we already assume the path provided
            else:
                
                oms.MzMLFile().load(self.input_path, self.input_map)
                self._log(
                    'Reading centroided data from `%s`.' % self.input_path
                )
    
    
    def find_features(self):
        
        self.do_mass_trace_detection()
        self.do_elution_peak_detection()
        self.do_feature_finding()
        self._log('Feature finding ready.')
    
    
    def do_mass_trace_detection(self):
        
        self._log('Performing mass trace detection.')
        self.mass_traces = []
        self.mass_trace_detection = MassTraceDetection(
            input_map = self.input_map,
            mass_traces = self.mass_traces,
            **self.mass_trace_detection_param
        )
        self.mass_trace_detection.main()
    
    
    def do_elution_peak_detection(self):
        
        self._log('Performing elution peak detection.')
        self.mass_traces_split = []
        self.elution_peak_detection = ElutionPeakDetection(
            mass_traces = self.mass_traces,
            mass_traces_split = self.mass_traces_split,
            **self.elution_peak_detection_param
        )
        self.elution_peak_detection.main()
    
    
    def do_feature_finding(self):
        
        self._log('Creating features by `FeatureFindingMetabo`.')
        self.output_map = oms.FeatureMap()
        self.mass_traces_filtered = []
        self.feature_finding_metabo = self.openms_obj
        self.adjust_featurefindermetabo_param()
        self.feature_finding_metabo.run(
            self.mass_traces_split,
            self.output_map,
            self.mass_traces_filtered,
        )
    
    
    def adjust_featurefindermetabo_param(self):
        
        self.param.pop(b'noise_threshold_int')
        self.param.pop(b'chrom_peak_snr')
        self.set_param()
    
    
    def write(self):
        """
        Writes the feature map into a ``featureXML`` file.
        """
        
        oms.FeatureXMLFile().store(self.output_path, self.output_map)
        self._log('Features have been written to `%s`.' % self.output_path)
    
    
    @staticmethod
    def _combine_param(common_param, param):
        
        combined_param = copy.deepcopy(common_param) or {}
        combined_param.update(param or {})
        
        return combined_param


class MapAlignmentAlgorithmPoseClustering(OpenmsMethodWrapper):
    """
    Applies a transformation on a feature map along its RT and m/z dimensions
    in order to remove noise and align traces from identical metabolites.
    
    Wrapper around ``pyopenms.MapAlignmentAlgorithmPoseClustering``.
    
    A map alignment algorithm based on pose clustering.

    Pose clustering analyzes pair distances to find the most probable
    transformation of retention times.
    The algorithm chooses the x most intensive peaks/features per map.
    This is modeled via the parameter ``max_num_peaks_considered``,
    which in turn influences the runtime and stability of the results.
    Bigger values prolong computation, smaller values might lead to no or
    unstable trafos. Set to -1 to use all features (might take very long for
    large maps).
    
    Parameters
    ----------
    input_file : str
        Path to a featureXML file.
    input_map : pyopenms.FeatureMap
        A``pyopenms.FeatureMap`` object. If provided, the data will be used
        from this instead of reading from ``input_file``. In this case
        ``input_file`` can be even ``None`` if ``output_file`` provided.
        Otherwise ``input_file`` will be used to create the ``output_file``
        path and name.
    output_file : str
        The output featureCML file where the aligned features will be saved
        If not provided it will be set according to the current settings.
    reference_file : str
        Path to the reference featureXML. The features will be aligned to
        this map.
    reference_map : pyopenms.FeatureMap
        A ``FeatureMap`` object instead of the file above. This way the
        reference features do not need to be read again.
    **kwargs
        Passed directly to ``pyopenms.MapAlignmentAlgorithmPoseClustering``.
    """
    
    
    def __init__(
            self,
            input_path = None,
            input_map = None,
            output_path = None,
            reference_path = None,
            reference_map = None,
            sample_id = None,
            sample_id_method = None,
            **kwargs
        ):
        
        OpenmsMethodWrapper.__init__(
            self,
            method = oms.MapAlignmentAlgorithmPoseClustering,
            input_path = input_path,
            output_path = output_path,
            sample_id = sample_id,
            sample_id_method = sample_id_method,
            name = 'map_aligner',
            **kwargs,
        )
        
        self.input_map = input_map
        self.reference_path = reference_path
        self.reference_map = reference_map
    
    
    def run(self):
        
        self.read_input()
        self.read_reference()
        self.align()
        self.write()
    
    
    def read_input(self):
        
        self._log('Obtaining input feature map.')
        
        self.input_map = self.read(
            path = self.input_path,
            mapobject = self.input_map,
        )
    
    
    def read_reference(self):
        
        if self.reference_path:
            
            self._log(
                'Reference feature map is from `%s`.' % self.reference_path
            )
        
        if not isinstance(self.reference_map, oms.FeatureMap):
            
            self._log('Obtaining reference feature map.')
            
            self.reference_map = self.read(
                path = self.reference_path,
                mapobject = self.reference_map,
            )
            
        else:
            
            self._log(
                'Reference map: using the provided '
                '`pyopenms.FeatureMap` object.'
            )
    
    
    def read(self, path, mapobject):
        
        if not isinstance(mapobject, oms.FeatureMap):
            
            mapobject = oms.FeatureMap()
            oms.FeatureXMLFile().load(path, mapobject)
            self._log('Feature map has been read from `%s`.' % path)
            
        else:
            
            self._log('Using the provided `pyopenms.FeatureMap` object.')
        
        return mapobject
    
    
    def align(self):
        
        self._log('Aligning the input map to the reference map.')
        # set reference_map file
        self.openms_obj.setReference(self.reference_map)
        # create object for the computed transformation
        self.trafo = oms.TransformationDescription()
        # do the alignment
        self.openms_obj.align(self.input_map, self.trafo)
        mat = oms.MapAlignmentTransformer()
        mat.transformRetentionTimes(self.input_map, self.trafo, False)
        dataproc = self.input_map.getDataProcessing()
        proc = oms.DataProcessing()
        proc.setProcessingActions({oms.ProcessingAction.ALIGNMENT})
        sw = proc.getSoftware()
        sw.setName('OpenMS')
        sw.setVersion(oms.VersionInfo.getVersion())
        proc.setSoftware(sw)
        proc.setCompletionTime(oms.DateTime.now())
        dataproc.append(proc)
        self.input_map.setDataProcessing(dataproc)
    
    
    def write(self):
        
        self.output_map = self.input_map
        # here we store the input_map as the transformation has been
        # applied to this map, the output_map is empty
        output_featurexml = oms.FeatureXMLFile()
        output_featurexml.store(self.output_path, self.input_map)
        self._log(
            'Aligned features have been written to `%s`.' % self.output_path
        )


class MgfExport(MethodPathHandler):
    """
    Exports the MS2 spectra in mascot generic format (MGF).
    
    Parameters
    ----------
    input_file : str
        Path to the mz
    """
    
    def __init__(
        self,
        input_path = None,
        input_obj = None,
        output_path = None,
        sample_id = None,
        sample_id_method = None,
    ):
        
        MethodPathHandler.__init__(
            self,
            input_path = input_path,
            input_obj = input_obj,
            output_path = output_path,
            name = 'mgf_export',
            method_key = 'mgf_export',
            sample_id = sample_id,
            sample_id_method = sample_id_method,
        )
    
    
    def main(self):
        
        self.set_paths()
        self.read()
        self.write()
    
    
    def read(self):
        
        if not isinstance(self.input_obj, oms.MSExperiment):
            
            mzml_file = oms.MzMLFile()
            self.input_obj = oms.MSExperiment()
            mzml_file.load(self.input_path, self.input_obj)
            self._log('Reading MS2 spectra from `%s`.' % self.input_path)
            
        else:
            
            self._log(
                'Using MS2 spectra from `pyopenms.MSExperiment` object.'
            )
    
    
    def write(self):
        
        with open(self.output_path, 'w') as fp:
            
            self._log('Exporting MS2 spectra to `%s`.' % self.output_path)
            # I commented this out at the moment,
            # I am not sure we need it:
            # Create header
            #outfile.write("COM=Testfile\n")
            #outfile.write("ITOL=1\n")
            #outfile.write("ITOLU=Da\n")
            #outfile.write("CLE=Trypsin\n")
            #outfile.write("CHARGE=1,2,3\n")
            
            # Iterate through all spectra,
            # skip all MS1 spectra and then write mgf format
            nr_ms2_spectra = 0
            
            for spectrum in self.input_obj:
                
                if spectrum.getMSLevel() == 1:
                    continue
                
                nr_ms2_spectra += 1
                _ = fp.write("BEGIN IONS\n")
                _ = fp.write("TITLE=%s\n" % spectrum.getNativeID().decode())
                _ = fp.write("RTINSECONDS=%.09f\n" % spectrum.getRT())
                self.spectrum = spectrum
                
                try:
                    _ = fp.write(
                        "PEPMASS=%.018f %.09f\n" % (
                            spectrum.getPrecursors()[0].getMZ(),
                            spectrum.getPrecursors()[0].getIntensity(),
                        )
                    )
                    ch = spectrum.getPrecursors()[0].getCharge()
                    
                    if ch > 0:
                        _ = fp.write("CHARGE=%u%s\n" % (
                                abs(ch),
                                '-' if ch < 0 else ''
                            )
                        )
                    
                except IndexError:
                    
                    _ = fp.write("PEPMASS=unknown\n")
                
                for peak in spectrum:
                    _ = fp.write(
                        "%s %s\n" % (peak.getMZ(), peak.getIntensity())
                    )
                
                _ = fp.write("END IONS\n")
            
            if nr_ms2_spectra == 0:
                
                self._log(
                    'Could not find any MS2 spectra in the input, '
                    'the output MGF file is empty!',
                    -1,
                )
                
            else:
                
                self._log(
                    '%u spectra have been written to `%s`.' % (
                        nr_ms2_spectra,
                        self.output_path,
                    )
                )


class FeatureGroupingAlgorithmQT(OpenmsMethodWrapper):
    """
    Wrapper around ``pyopenms.FeatureGroupingAlgorithmQT``.
    """
    
    
    def __init__(
            self,
            input_path = None,
            input_obj = None,
            output_path = None,
            sample_id = None,
            sample_id_method = None,
            **kwargs,
        ):
        
        OpenmsMethodWrapper.__init__(
            self,
            method = oms.PeakPickerHiRes,
            input_path = input_path,
            input_obj = input_obj,
            output_path = output_path,
            name = 'peak_picker',
            sample_id = sample_id,
            sample_id_method = sample_id_method,
            **kwargs
        )


    def main(
            self,
            input_path = None,
            output_path = None,
            keep_subelements = True,
            **params
        ):
        
        #------------------------------------------
        #   service utils;
        #------------------------------------------
        def addDataProcessing(obj, params, action):
            if isinstance(obj, oms.MSExperiment):
                result = oms.MSExperiment()
                for spec in obj:
                    spec = _addDataProcessing(spec, params, action)
                    result.addSpectrum(spec)
            else:
                result = _addDataProcessing(obj, params, action)
            return result
        
        
        def _addDataProcessing(item, params, action):
            dp = item.getDataProcessing()
            p = oms.DataProcessing()
            p.setProcessingActions(set([action]))
            sw = p.getSoftware()
            sw.setName(os.path.basename(sys.argv[0]))
            sw.setVersion(oms.VersionInfo.getVersion())
            p.setSoftware(sw)
            p.setCompletionTime(oms.DateTime.now())

            for k, v in params.items():
                p.setMetaValue("parameter: "+k, v)

            dp.append(p)
            item.setDataProcessing(dp)
            return item
        
        
        def writeParamsIfRequested(args, params):
            if args.write_dict_ini or args.write_ini:
                if args.write_dict_ini:
                    with open(args.write_dict_ini, "w") as fp:
                        pprint.pprint(params.asDict(), stream=fp)
                if args.write_ini:
                    fh = oms.ParamXMLFile()
                    fh.store(args.write_ini, params)
                return True
            return False
        
        
        def updateDefaults(args, defaults):
            if args.ini:
                param = oms.Param()
                fh = oms.ParamXMLFile()
                fh.load(args.ini, param)
                defaults.update(param)
            elif args.dict_ini:
                with open(args.dict_ini, "r") as fp:
                    try:
                        dd = eval(fp.read())
                    except:
                        raise Exception("could not parse %s" % args.dict_ini)
                defaults.update(dd)
        
        
        def link(input_path = None, output_path = None, keep_subelements = True, **params):
            if not (input_path and output_path) :
                raise Exception("source file pattern have not be empty!")
            
            from  collections import Counter

            in_files = get_list_full_names( input_path )
            out_file = output_path
            
            in_types = set(oms.FileHandler.getType(in_) for in_ in in_files)
            
            link_features = True
            """
            if in_types == set((oms.Type.CONSENSUSXML,)):
                link_features = False
            elif in_types == set((oms.Type.FEATUREXML,)):
                link_features = True
            else:
                raise Exception("different kinds of input files")
            """
            #algorithm_parameters = params.copy("algorithm:", True)
            #algorithm_parameters = {  : }
            algorithm = oms.FeatureGroupingAlgorithmQT()
            p = algorithm.getDefaults()
            p.setValue("algorithm", True)
            #algorithm.setParameters(algorithm_parameters)
            algorithm.setParameters( p )

            out_map = oms.ConsensusMap()
            fds = out_map.getColumnHeaders()
            if link_features:
                f = oms.FeatureXMLFile()
                maps = []
                for i, in_file in enumerate(in_files):
                    map_ = oms.FeatureMap()
                    f.load(in_file, map_)

                    # set filedescriptions
                    fd = fds.get(i, oms.ColumnHeader())
                    fd.filename = in_file
                    fd.size = map_.size()
                    fd.unique_id = map_.getUniqueId()
                    fds[i] = fd
                    maps.append(map_)
                out_map.setColumnHeaders(fds)
                algorithm.group(maps, out_map)
            else:
                f = oms.ConsensusXMLFile()
                maps = []
                for i, in_file in enumerate(in_files):
                    map_ = oms.ConsensusMap()
                    f.load(in_file, map_)
                    maps.append(map_)
                algorithm.group(maps, out_map)

                if not keep_subelements:
                    for i in range(len(in_files)):
                        # set filedescriptions
                        fd = fds.get(i, oms.ColumnHeader())
                        fd.filename = in_files[i]
                        fd.size = maps[i].size()
                        fd.unique_id = maps[i].getUniqueId()
                        fds[i] = fd
                    out_map.setColumnHeaders(fds)
                else:
                    algorithm.transferSubelements(maps, out_map)

            out_map.setUniqueIds()
            addDataProcessing(out_map, params, oms.ProcessingAction.FEATURE_GROUPING)

            oms.ConsensusXMLFile().store(out_file, out_map)
            """
            sizes = []
            for feat in out_map:
                sizes.append(feat.size())

            c = Counter(sizes)
            print "Number of consensus features:"
            for size, count in c.most_common():
                print "   of size %2d : %6d" % (size, count)
            print "        total : %6d" % out_map.size()
        #   service utils;
            """
        # --- end of link ---
        
        link(input_path = input_path, output_path = output_path, keep_subelements = keep_subelements, **params)
    # --- end of main2 ---


class MSPreprocess(PathHandlerBase):
    """
    
    Workflow covering all MS preprocessing steps: peak picking, feature
    construction, map alignment and feature grouping. Also exports MS2
    spectra in MGF format.
    
    Parameters
    ----------
    force : set
        A set of steps which must be done no matter if the output of a later
        stage is available. E.g. ``{'peak_picking'}``.
    stop : str
        A step where the workflow has to stop. E.g. if it's
        ``feature_finding``, the last step executed will be
        ``feature_finding``.
    """
    
    _stages = (
        'profile',
        'centroided',
        'features',
        'features_aligned',
        'features_grouped',
    )
    _steps = (
        ('profile', 'centroided', 'peak_picking'),
        ('centroided', 'features', 'feature_finding'),
        ('features', 'features_aligned', 'map_alignment'),
        ('features_aligned', 'features_grouped', 'feature_grouping'),
    )
    _step_data = dict(
        (
            method,
            (source, target)
        )
        for source, target, method in _steps
    )
    _stages_by_source = dict(
        (
            source,
            method
        )
        for source, target, method in _steps
    )
    _stages_by_source['features_grouped'] = None
    _classes = {
        'peak_picking': PeakPickerHiRes,
        'feature_finding': FeatureFindingMetabo,
        'map_alignment': MapAlignmentAlgorithmPoseClustering,
        'feature_grouping': FeatureGroupingAlgorithmQT,
    }
    _input_extensions = {
        'profile': 'mzML',
        'centroided': 'mzML',
        'features': 'featureXML',
        'features_aligned': 'featureXML',
        'features_grouped': 'consensusXML',
    }
    
    def __init__(
        self,
        # stages
        profile = None,
        centroided = None,
        features = None,
        features_aligned = None,
        features_grouped = None,
        # further parameters for input/output
        output_dir = None,
        output_path = None,
        sample_ids = None,
        sample_id_method = None,
        input_path = None,
        input_obj = None,
        input_filter = None,
        # parameters for each step
        peak_picking_param = None,
        mass_trace_detection_param = None,
        elution_peak_detection_param = None,
        feature_finding_metabo_param = None,
        feature_finding_common_param = None,
        map_alignment_param = None,
        reference_sample = None,
        feature_grouping_param = None,
        # workflow parameters
        force = None,
        stop = None,
        mgf_export = True,
        # nothing
        **kwargs
    ):
        
        session.Logger.__init__(self, name = 'msproc')
        
        self.method_key = 'ms_preprocess'
        self.peak_picking_param = peak_picking_param or {}
        self.mass_trace_detection_param = mass_trace_detection_param or {}
        self.elution_peak_detection_param = elution_peak_detection_param or {}
        self.feature_finding_metabo_param = feature_finding_metabo_param or {}
        self.feature_finding_common_param = feature_finding_common_param or {}
        self.map_alignment_param = map_alignment_param or {}
        self.feature_grouping_param = feature_grouping_param or {}
        
        self.sample_ids = sample_ids
        self.force = force or set()
        self.profile = profile
        self.centroided = centroided
        self.features = features
        self.features_aligned = features_aligned
        self.features_grouped = features_grouped
        
        self._input_path = input_path
        self._input_obj = input_obj
        self.input_filter = input_filter
        self.sample_id_method = sample_id_method
        self.start_from = 'profile'
        
        self._set_inputs()
        
        PathHandlerBase.__init__(
            self,
            input_path = self._input_path,
            input_obj = self._input_obj,
            input_filter = self.input_filter,
            input_ext = self._input_extensions[self.start_from],
            sample_id_method = self.sample_id_method,
            sample_id = self.sample_ids,
            method_key = self.start_from,
        )
        PathHandlerBase.set_paths(self)
        
        _ = self._ensure_tuple('_input_path', common.basestring)
        _ = self._ensure_tuple('_input_obj', OPENMS_OBJ_TYPES)
        
        self.multi_file_input = True
        self._set_sample_id()
        
        self._set_first_input()
        self.mgf_export = mgf_export
        self.reference_sample = reference_sample
    
    
    def reload(self):
        modname = self.__class__.__module__
        mod = __import__(modname, fromlist = [modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)
    
    
    def main(self):
        
        self.workflow()
    
    
    def workflow(self):
        
        PathHandlerBase.set_paths(self)
        self._set_sample_ids()
        
        for source, result, method in self._steps:
            
            if not getattr(self, result) or method in self.force:
                
                self._log('Next step in workflow: `%s`.' % method)
                getattr(self, method)()
            
            if method == self.stop:
                
                self._log('Stopping after step `%s`.' % method)
                break
    
    
    def _ensure_tuple(self, attr, types):
        """
        Makes sure an attribute is tuple (or list or array) even if a
        single element one.
        Returns a boolean value which confirms if the attribute is indeed
        is a tuple like object.
        """
        
        if isinstance(getattr(self, attr), types):
            
            setattr(self, attr, (getattr(self, attr),))
        
        return isinstance(getattr(self, attr), (tuple, list, np.ndarray))
    
    
    @staticmethod
    def _check_type(obj, types):
        
        return (
            isinstance(obj, types) or (
                isinstance(obj, (tuple, list, np.ndarray)) and
                all(isinstance(o, types) for o in obj)
            )
        )
    
    
    def _set_inputs(self):
        
        _input_types = (common.basestring,) + OPENMS_OBJ_TYPES
        
        # if we force to re-do from the first possible step we iterate
        # from the beginning, else we go backwards to find the last
        # available stage
        
        if all(getattr(self, stage) is None for stage in self._stages):
            
            setattr(self, self._stages[0], True)
        
        stages = reversed(self._stages)
        stages_covered = set()
        
        self._inputs = ()
        
        for stage in stages:
            
            stages_covered.add(self._stages_by_source[stage])
            
            # this will be the first stage if...
            if (
                # we have this stage set to True or input data provided
                getattr(self, stage) and
                # and executing all steps not forced
                self.force != True and
                # and all forced steps already covered
                not self.force - stages_covered
            ):
                
                # it's a boolean True, meaning that we start from this
                # stage and get the inputs from `input_path` or `input_obj`
                # arguments
                if getattr(self, stage) == True:
                    
                    if self._input_obj is not None:
                        
                        setattr(self, stage, self._input_obj)
                        
                    elif self._input_path is not None:
                        
                        setattr(self, stage, self._input_path)
                
                self.start_from = stage
                first_stage = getattr(self, stage)
                
                # looks like a single path or object provided
                if self._check_type(stage, OPENMS_OBJ_TYPES):
                    
                    self._input_obj = first_stage
                
                elif self._check_type(stage, common.basestring):
                    
                    self._input_path = first_stage
                
                break
    
    
    def _set_first_input(self):
        
        _input = (
            self._input_obj
                if self._input_obj is not None else
            self._input_path
        )
        setattr(self, self.start_from, _input)
    
    
    def _set_sample_ids(self):
        
        if not self.sample_ids:
            
            self.sample_ids = [None] * len(self._inputs)
    
    
    def _step_base(self, method):
        
        source_name, target_name = self._step_data[method]
        source = getattr(self, source_name)
        source_path_attrs = '%s_paths' % source_name
        source_paths = (
            getattr(self, source_path_attrs)
                if hasattr(self, source_path_attrs) else
            [None] * len(source)
        )
        target = []
        target_paths = []
        _class = self._classes[method]
        param = getattr(self, '%s_param' % method)
        
        sample_ids_in = self.sample_id or [None] * len(source)
        sample_ids_out = []
        
        for resource, source_path, sample_id in zip(
                source,
                source_paths,
                sample_ids_in,
            ):
            
            input_arg = (
                'input_obj'
                    if isinstance(resource, OPENMS_OBJ_TYPES) else
                'input_path'
            )
            param = copy.deepcopy(param)
            param[input_arg] = resource
            param['sample_id'] = sample_id
            
            if input_arg != 'input_path':
                
                param['input_path'] = source_path
            
            worker = _class(**param)
            worker.main()
            target.append(worker.output_map)
            target_paths.append(worker.output_path)
            sample_ids_out.append(worker.sample_id)
            
            if self.mgf_export and method == 'peak_picking':
                
                self.export_mgf(
                    input_obj = worker.output_map,
                    input_path = worker.output_path,
                )
        
        self.result = target
        self.result_paths = target_paths
        setattr(self, target_name, target)
        setattr(self, '%s_paths' % target_name, target_paths)
        # if we already had sample IDs in the previous step
        # they will just pass through
        self.sample_ids = sample_ids_out
    
    
    def peak_picking(self):
        
        self._step_base(method = 'peak_picking')
    
    
    def export_mgf(self, **kwargs):
        
        mgf_exporter = MgfExport(**kwargs)
        mgf_exporter.main()
    
    
    def feature_finding(self):
        
        self._setup_feature_finding_param()
        
        self._step_base(method = 'feature_finding')
    
    
    def _setup_feature_finding_param(self):
        
        self.feature_finding_param = copy.deepcopy(
            self.feature_finding_common_param
        )
        self.feature_finding_param['mass_trace_detection_param'] = (
            self.mass_trace_detection_param
        )
        self.feature_finding_param['elution_peak_detection_param'] = (
            self.elution_peak_detection_param
        )
        self.feature_finding_param['feature_finding_metabo_param'] = (
            self.feature_finding_metabo_param
        )
    
    
    def map_alignment(self):
        
        # attempting to select if any sample requested
        if self.reference_sample is not None:
            
            self._map_alignment_manually_select_reference()
        
        # falling back to selecting the largest map
        if (
            'reference_map' not in self.map_alignment_param and
            'reference_path' not in self.map_alignment_param
        ):
            self._map_alignment_largest_map_as_reference()
        
        self._step_base(method = 'map_alignment')
    
    
    def _map_alignment_largest_map_as_reference(self):
        
        sizes = []
        
        for i, femap in enumerate(self.features):
            
            if isinstance(femap, common.basestring) and os.path.exists(femap):
                
                fh = oms.FeatureXMLFile()
                sizes.append(fh.loadSize(femap))
                
            else:
                
                sizes.append(femap.size())
        
        iref = self.sample_id[sizes.index(max(sizes))]
        refmap = self.features[iref]
        
        self._log(
            'Selected the largest sample as reference for the alignment: '
            'sample `%s` with %u features.' % (
                str(self.sample_id[iref]),
                sizes[iref],
            )
        )
        
        if isinstance(refmap, oms.FeatureMap):
            
            self.map_alignment_param['reference_map'] = refmap
            
        else:
            
            self.map_alignment_param['reference_path'] = refmap
    
    
    def _map_alignment_manually_select_reference(self, sample_id = None):
        
        reference_sample = sample_id or self.reference_sample
        
        if reference_sample:
            
            if reference_sample not in self.sample_id:
                
                self._log(
                    'Could not select reference sample `%s`. '
                    'Sample ID does not exist.' % str(reference_sample)
                )
                
            else:
                
                self._log(
                    'Manually selected reference sample: `%s`.' % (
                        self.reference_sample
                    )
                )
                idx = self.sample_id.index(self.reference_sample)
                self.map_alignment_param['reference_map'] = self.result[idx]
                self.map_alignment_param['reference_path'] = (
                    self.result_paths[idx]
                )
    
    
    def feature_grouping(self):
        
        self._step_base(method = 'feature_grouping')
    
    
    def get_sampleset(self):
        
        return {}
