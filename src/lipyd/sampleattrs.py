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

from past.builtins import xrange, range

import numpy as np

import lipyd.common as common


class SampleAttrs(object):
    
    def __init__(self, sample_id_method = None, attrs = None):
        """
        Represents the ID and attributes of a sample.
        
        :param object,callable sample_id_method:
            Either an object (string or tuple) or a method which takes
            a dict of attributes as its argument and returns a unique
            identifier for the sample.
        :param dict attrs:
            A dictionary of sample attributes.
        """
        
        self.attrs = attrs or {}
        self._sample_id_method = sample_id_method
        self._set_sample_id()
    
    def _set_sample_id(self):
        
        if self._sample_id_method is None:
            
            # first if it's None we call the deafult method to
            # create sample ID from the sample attributes
            return self._default_sample_id_method()
            
        elif callable(self._sample_id_method):
            
            # if a custom method has been provided we use
            # that instead
            return self._sample_id_method(self.attrs)
            
        else:
            
            # if it's not callable but any other kind of object
            # then we assume the sample ID is explicitely given
            return self._sample_id_method
    
    def _default_sample_id_method(self):
        
        if (
            self.attrs is not None and
            'label' in self.attrs and
            'fraction' in self.attrs['label']
        ):
            
            return self.attrs['label']['fraction']
        
        return common.random_string()


class SampleSetAttrs(object):
    
    def __init__(
            self,
            sample_ids_method = None,
            sample_attrs = None,
            sample_id_processor = None,
            length = None,
        ):
        
        # custom method to process sample IDs
        self._sample_id_processor = sample_id_processor
        
        if isinstance(sample_ids_method, (list, np.ndarray)):
            
            # sample IDs provided for each sample
            sample_ids = [
                self._get_sample_id(sample_id)
                for sample_id in sample_ids
            ]
            
            length = len(sample_ids)
            
        elif sample_attrs is not None:
            
            # if sample attributes provided
            length = len(sample_attrs)
        
        # now it is either None or a method
        sample_ids = [sample_ids_method] * length
        
        if length is None:
            
            raise RuntimeError(
                'SampleSetAttrs: number of samples not provided.'
            )
        
        self.attrs = np.array([
            SampleAttrs(
                sample_id_method = sample_ids[i],
                attrs = sample_attrs[i] if sample_attrs else {},
            )
            for i in xrange(length)
        ])
        
        self._set_sample_ids()
    
    def _get_sample_id(self, sample_id):
        
        if callable(self._sample_id_processor):
            
            return self._sample_id_processor(sample_id)
        
        return sample_id
    
    def _set_sample_ids(self):
        # Note: this is called by Sample.__init__()
        
        self.sample_index_to_id = []
        
        for attr in self.attrs:
            
            self.sample_index_to_id.append(attr.sample_id)
        
        self._update_id_to_index()
    
    def _update_id_to_index(self):
        
        self.sample_id_to_index = dict(
            reversed(i)
            for i in enumerate(self.sample_index_to_id)
        )
    
    def get_sample_id(self, i):
        """
        Returns the ID of a sample by its index.
        """
        
        return self.sample_index_to_id[i]
    
    def get_sample_index(self, sample_id):
        """
        Returns the index of a sample by its ID.
        """
        
        return self.sample_id_to_index[sample_id]
    
    def sort_by_index(self, idx):
        
        idx = np.array(idx)
        
        self.attrs = self.attrs[idx]
        
        self._set_sample_ids()
    
    def sort_by_sample_id(self, sample_ids, process = False):
        """
        Sets the order of the sample attributes according to the list of
        sample IDs provided.
        
        :param list sample_ids:
            A list of sample IDs, e.g. ``[('A', 1), ('A', 2), ...]``.
        :param bool process:
            Process ``sample_ids`` by the ``sample_id_processor`` method.
        """
        
        idx = np.array([
            self.sample_id_to_index[
                sample_id if not process else self._get_sample_id(sample_id)
            ]
            for sample_id in sample_ids
        ])
        
        self.sort_by_index(idx)
