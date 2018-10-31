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

from future.utils import iteritems

import numpy as np

import lipyd.sampleattrs as sampleattrs


class SampleData(sampleattrs.SampleSorter):
    
    def __init__(
            self,
            samples = None,
            sample_ids = None,
            sample_data = None,
            sample_id_processor = None,
            **kwargs,
        ):
        """
        Represents data about a series of samples. Samples are LC MS/MS runs.
        Data might be a binary, qualitative or quantitative attribute about
        these samples. E.g. disease status of a patient, time of sampling,
        quantity of a protein, etc. The data might have more than one
        dimensions but the first axis is always considered to be the sample
        identity and the number of samples must agree with the number of
        samples in the sampleset. If no ``sample_ids`` provided the order of
        sample data assumed to be the same as the order of samples in the
        sampleset. Otherwise it will be ordered according to the labels.
        
        :param 
        :param list,numpy.array **kwargs:
            Data associated to samples. First dimension of each array must
            agree with the samples in the sampleset. Multiple variables
            might be provided each will be set as an attribute of the object.
        :param lipyd.samples.SampleSet samples:
            A ``SampleSet`` object.
        :param list sample_ids:
            A list of sample IDs. Will be used to reorder the sample data
            in order to have the same ordering as the samples in sampleset.
        :param callable sample_id_processor:
            A method to process elements in ``sample_ids``.
        """
        
        if (not sample_ids and not samples) or not hasattr(samples, 'attrs'):
            
            raise RuntimeError(
                'SampleData: either `samples` or '
                '`sample_ids` must be provided.'
            )
        
        sample_ids = sample_ids or samples.attrs.sample_index_to_id
        
        self.numof_samples = len(sample_ids)
        
        if isinstance(samples, sampleattrs.SampleSorter):
            
            sample_data = sample_data or []
            sample_data.append(samples)
        
        self.var = set()
        
        for attr, data in iteritems(kwargs):
            
            self._add_var(data, attr)
        
        sampleattrs.SampleSorter.__init__(
            self,
            sample_data = sample_data,
            sample_ids = sample_ids,
            sample_axis = 0,
            sample_id_processor = sample_id_processor,
        )
    
    def _add_var(self, data, attr):
        """
        Adds a new variable to the data handler object.
        The first dimension of the array must agree with
        the number of samples.
        """
        
        if isinstance(data, list):
            
            data = np.array(data)
        
        if data.shape[0] != self.numof_samples:
            
            raise RuntimeError(
                'SampleData: first dimension of each array must agree with '
                'the number of samples. `%s` has length %u while numof '
                'samples is %u.' % (
                    attr, data.shape[0], self.numof_samples,
                )
            )
        
        setattr(self, attr, data)
        self.var.add(attr)
    
    @staticmethod
    def _bool_array(
            selection,
            samples = None,
            sample_ids = None,
            sample_id_processor = None,
        ):
        """
        This method helps to set up a sample selection but placed here as
        staticmethod to make it usable other places.
        """
        
        if not isinstance(selection[0], bool):
            
            if sample_id_processor:
                
                selection = [sample_id_processor(s) for s in selection]
            
            if not sample_ids:
                
                if isinstance(samples, sampleattrs.SampleSorter):
                    
                    first = samples
                    
                else:
                    
                    first = sample_data[0]
                
                sample_ids = first.attrs.sample_index_to_id
            
            selection = np.array([s in selection for s in sample_ids])
        
        return selection


class SampleSelection(SampleData):
    
    def __init__(
            self,
            selection,
            samples = None,
            sample_ids = None,
            sample_data = None,
            sample_id_processor = None,
            **kwargs,
        ):
        """
        Represents a binary selection of samples. You can provide either a
        boolean array or a list of sample IDs to be selected.
        """
        
        sample_id_processor = sample_id_processor or (lambda x: x)
        
        selection = self._bool_array(
            selection = selection,
            samples = samples,
            sample_ids = sample_ids,
            sample_id_processor = sample_id_processor,
        )
        
        SampleData.__init__(
            self,
            samples = samples,
            sample_ids = sample_ids,
            sample_data = sample_data,
            selection = selection,
            sample_id_processor = sample_id_processor,
        )
    
    def sample_ids_selected(self):
        """
        Returns a list IDs of selected samples.
        """
        
        return [
            sample_id
            for i, sample_id in self.attrs.sample_index_to_id
            if self.selection[i]
        ]


class FeatureAnalyzer(object):
    
    def __init__(self, samples):
        
        self.samples = samples


class ProfileFeatureAnalyzer(FeatureAnalyzer):
    
    def __init__(self, samples):
        
        FeatureAnalyzer.__init__(self, samples = samples)
