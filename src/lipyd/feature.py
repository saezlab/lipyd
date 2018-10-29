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

import lipyd.sampleattrs as sampleattrs


class SampleData(sampleattrs.SampleSorter):
    
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
        
        SampleSorter.__init__(
            self,
            sample_data = sample_data,
            sample_axis = 0,
        )
    
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
        
        


class SampleSelection(SampleData):
    
    def __init__(self, selection, labels = None):
        
        pass


class FeatureAnalyzer(object):
    
    def __init__(self, samples):
        
        self.samples = samples


class ProfileFeatureAnalyzer(FeatureAnalyzer):
    
    def __init__(self, samples):
        
        FeatureAnalyzer.__init__(self, samples = samples)
