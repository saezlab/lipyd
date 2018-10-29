#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `lipyd` python module
#
#  Copyright (c) 2015-2018 - EMBL
#
#  File author(s): Dénes Türei (turei.denes@gmail.com)
#
#  Distributed under the GNU GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://www.ebi.ac.uk/~denes
#

import numpy as np


class SampleSorter(object):
    
    def __init__(self, sample_data = None, sample_axis = 0):
        """
        Keeps the order of samples synchronized between multiple objects.
        These objects represent sets of the same samples such as
        ``sample.SampleSet`` or ``feature.SampleData``.
        
        :param list,set sample_data:
            Other ``sample.SampleSet`` or ``feature.SampleData`` derived
            objects that should keep the same order of samples.
        :param int sample_axis:
            Which axis in the arrays corresponds to the samples.
            In ``sample.SampleSet`` objects this is axis 1 as axis 0
            corresponds to the features. In ``feature.SampleData`` derived
            objects this is axis 0.
        """
        
        self.sample_data = {}
        self._sample_axis = sample_axis
        
        if sample_data is None:
            
            sample_data = []
        
        if not isinstance(sample_data, (list, set)):
            
            sample_data = [sample_data]
        
        for s in sample_data:
            
            self.register(s)
    
    def register(self, s):
        
        self.sample_data[id(s)] = s
    
    def order_samples(self, idx, done = None):
        """
        Changes the ordering of the samples.
        """
        
        done = set() if done is None else done
        
        if id(self) in done:
            
            return
        
        numof_samples = self.numof_samples
        
        if len(idx) != numof_samples:
            
            raise RuntimeError(
                'Invalid index length: %u while number of samples is %u.' % (
                    len(idx), numof_samples
                )
            )
        
        for var in self.var:
            
            arr = getattr(self, var)
            
            if len(arr.shape) < 1 or arr.shape[1] != numof_samples:
                
                continue
            
            setattr(self, var, np.take(arr, idx, axis = 1))
        
        done.add(id(self))
        
        for sd_id, sd in iteritems(self.sample_data):
            
            if sd_id not in done:
                
                sd.order_samples(idx = idx, done = done)


class SampleData(SampleSorter):
    
    def __init__(
            self,
            data,
            samples = None,
            labels = None,
            sample_data = None,
        ):
        """
        Represents data about a series of samples. Samples are LC MS/MS runs.
        Data might be a binary, qualitative or quantitative attribute about
        these samples. E.g. disease status of a patient, time of sampling,
        quantity of a protein, etc. The data might have more than one
        dimensions but the first axis is always considered to be the sample
        identity and the number of samples must agree with the number of
        samples in the sampleset. If no ``labels`` provided the order of
        sample data assumed to be the same as the order of samples in the
        sampleset. Otherwise it will be ordered according to the labels.
        
        :param numpy.array data:
            Data associated to samples. First dimension must agree with the
            samples in the sampleset.
        :param lipyd.samples.SampleSet samples:
            A ``SampleSet`` object.
        :param list labels:
            A list of sample labels. Will be used to reorder the sample data
            in order to have the same ordering as the samples in sampleset.
        """
        
        self.data    = data
        self.labels  = labels
        self.samples = samples
        
        SampleSorter.__init__(self, sample_data = sample_data)
    
    def __len__(self):
        """
        Returns the number of samples.
        """
        
        return self.data.shape[0]
    
    def set_sampleset(self, sampleset):
        """
        Registers a ``SampleSet`` object.
        """
        
        if sampleset.numof_samples != len(self):
            
            raise RuntimeError(
                'First dimension of the data array '
                'must gree with the number of samples.'
            )
        
        


class SampleSelection(object):
    
    def __init__(self, selection, labels = None):
        
        pass


class FeatureAnalyzer(object):
    
    def __init__(self, samples):
        
        self.samples = samples


class ProfileFeatureAnalyzer(FeatureAnalyzer):
    
    def __init__(self, samples):
        
        FeatureAnalyzer.__init__(self, samples = samples)
