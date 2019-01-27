#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `lipyd` python module
#
#  Copyright (c) 2014-2019 - EMBL
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

import imp

import numpy as np


class Recalibration(object):
    
    
    def __init__(
            self,
            first = None,
            last = None,
            by_sample = None,
            **kwargs,
        ):
        """
        Provides recalibration data for ``lipyd.sample.Sample`` and
        ``lipyd.sample.SampleSet`` objects.
        
        :arg float,numpy.ndarray first:
            Error in ppm for the first sample. Either a single float or an
            array of floats with the same length as number of features in
            the sample or an array of 2 columns: first with m/z values,
            second with errors in ppm above the corresponding m/z value.
        :arg float,numpy.ndarray last:
            Error in ppm for the last sample. Same as ``first``. If not
            provided all samples will be corrected by ``first``.
        :arg dict by_sample:
            Errors for each sample. Keys are sample IDs, values are as
            detailed at argument ``first``: float or array with one or
            two columns. ``**kwargs`` handled the same way as ``by_sample``
            but the latter has priority.
        """
        
        self.first = first
        self.last = last
        self.by_sample = by_sample or kwargs
    
    
    def reload(self):
        """ """
        
        modname = self.__class__.__module__
        mod = __import__(modname, fromlist = [modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)
    
    
    def recalibrate(self, samples):
        """
        Performs correction on m/z values in a ``lipyd.sample.Sample`` or a
        ``lipyd.sample.SampleSet`` object.
        
        :arg Sample,SampleSet samples:
            A ``lipyd.sample.Sample`` or a ``lipyd.sample.SampleSet`` object.
        """
        
        # processing ppms
        self.samples = samples
        self.numof_samples = self.samples.numof_samples
        self.numof_features = len(self.samples)
        self.samples.sort_all('mzs')
        self.sample_id_proc = samples.attrs.proc
        self.process_ppms()
        self._by_feature_ppms()
        
        # keeping the original mzs
        self.samples._add_var(self.samples.mzs, 'mzs_original')
        self.samples._add_var(
            self.samples.mzs_by_sample,
            'mzs_by_sample_original'
        )
        
        # doing the recalibration
        self.samples.mzs_by_sample = self._recalibrate(
            mzs = self.samples.mzs_by_sample,
            ppms = self.ppms,
        )
        
        self.samples.mzs = self._recalibrate(
            mzs = self.samples.mzs,
            ppms = self.ppms_by_feature,
        )
    
    
    def _recalibrate(self, mzs, ppms):
        
        return mzs - mzs / 1e6 * ppms
    
    
    def process_ppms(self):
        """
        Arranges ppms in a way as described argument definitions for
        ``__init__``, taking into consideration the current samples object.
        """
        
        if self.first is not None:
            
            if isinstance(self.first, (float, int, np.float64)):
                
                if self.last is None:
                    
                    # we got one single ppm value for all samples
                    self.ppms = self.first
                    
                else:
                    
                    # we got 2 values for first and last sample, respectively
                    # extrapolate these over all samples
                    self.ppms = np.linspace(
                        self.first,
                        self.last,
                        self.numof_samples,
                    )
                
            elif isinstance(self.first, np.ndarray):
                
                # we got array of ppms for the first sample, this means
                # the ppms are not even across the m/z range
                # here we make sure we have an array of ppms with the
                # same length as number of features
                mzs_first = (
                    self.samples.mzs_by_sample[:,0]
                        if hasattr(self.samples, 'mzs_by_sample') else
                    self.samples.mzs
                )
                first = self.align_array(self.first, mzs_first)
                
                if not isinstance(self.last, np.ndarray):
                    
                    # ppm values provided only for the first sample
                    # we use this for all samples
                    self.ppms = first
                    
                else:
                    
                    # ppms for the last sample also provided
                    # we extrapolate for other samples in between
                    mzs_last = (
                        self.samples.mzs_by_sample[:,-1]
                            if hasattr(self.samples, 'mzs_by_sample') else
                        self.samples.mzs
                    )
                    last = self.align_array(self.last, mzs_last)
                    
                    ppms = []
                    
                    for f, l in zip(first, last):
                        
                        ppms.append(np.linspace(f, l, self.numof_samples))
                    
                    self.ppms = np.vstack(ppms)
            
        if (isinstance(self.by_sample, (list, np.ndarray)) and
            len(self.by_sample) == self.numof_samples
        ):
            
            # if it is a list or array we convert to dict so we can process
            # the same way as dicts
            self.by_sample = dict(
                zip(
                    self.samples.sample_index_to_id,
                    self.by_sample
                )
            )
        
        if self.by_sample is not None:
            
            self.ppms = []
            
            # we got ppms for each sample, iterate through the dict
            for sample_id, sample_ppm in iteritems(self.by_sample):
                
                if isinstance(sample_ppm, (float, int, np.float64)):
                    
                    self.ppms.append(sample_ppm)
                    
                else:
                    
                    sample_id = self.sample_id_proc(sample_id)
                    
                    sampl_mzs = mzs_last = (
                        self.samples.mzs_by_sample[
                            :,
                            self.samples.attrs.sample_id_to_index[sample_id]
                        ]
                        # this only to keep compatibility with Sample
                            if hasattr(self.samples, 'mzs_by_sample') else
                        self.samples.mzs
                    )
                    sample_ppm = self.align_array(sample_ppm, sample_mzs)
                    self.ppms.append(sample_ppm)
            
            if all(isinstance(ppms, np.ndarray) for ppms in self.ppms):
                
                # if above we produced ppms for each feature
                # we create an array with dimensions features x samples
                self.ppms = np.hstack(self.ppms)
                
            if isinstance(self.ppms, list):
                
                # if we have one number for each sample, simply
                # create a one dimensional array
                self.ppms = np.array(self.ppms)
    
    
    def _by_feature_ppms(self):
        
        if isinstance(self.ppms, (float, int, np.float64)):
            
            # if it's a single number
            self.ppms_by_feature = self.ppms
            
        elif isinstance(self.ppms, np.ndarray) and self.ppms.ndim == 1:
            
            # if it's an array with one number for each sample
            self.ppms_by_feature = np.median(self.ppms)
            
        else:
            
            # if it's an array with one line for each feature
            self.ppms_by_feature = np.median(self.ppms, 1)
    
    
    def align_array(self, array, mzs):
        
        if len(array) == self.numof_features:
            
            return array
            
        else:
            
            if ndim(array) < 2:
                
                raise(
                    ValueError,
                    'ppm arrays must be the same length as '
                    'number of features in the samples or must have '
                    '2 columns: first with m/z values second with ppms.'
                )
            
            # sorting by m/z
            array = array[array[:,0].argsort(),:]
            
            ppms = []
            
            ppm = array[0,1]
            i = 1
            
            for mz in mzs:
                
                if mz >= array[i,0]:
                    
                    ppm = array[i,1]
                    
                    if i < array.shape[0]:
                        
                        i += 1
                
                ppms.append(ppm)
            
            return np.array(ppms)
