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
        self.sample_id = self._get_sample_id()
    
    def _get_sample_id(self):
        
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
                for sample_id in sample_ids_method
            ]
            
            length = len(sample_ids)
            
        elif sample_attrs is not None:
            
            # if sample attributes provided
            length = len(sample_attrs)
            
        else:
            
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
    
    def __len__(self):
        
        return len(self.attrs)
    
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
    
    def argsort_by_sample_id(self, sample_ids, process = False):
        """
        Returns an index array which sorts the sample attributes according
        to the list of sample IDs provided.
        
        :param list sample_ids:
            A list of sample IDs, e.g. ``[('A', 1), ('A', 2), ...]``.
        :param bool process:
            Process ``sample_ids`` by the ``sample_id_processor`` method.
        """
        
        return np.array([
            self.sample_id_to_index[
                sample_id if not process else self._get_sample_id(sample_id)
            ]
            for sample_id in sample_ids
        ])
    
    def sort_by_sample_id(
            self,
            sample_ids,
            process = False,
            return_idx = False,
        ):
        """
        Sets the order of the sample attributes according to the list of
        sample IDs provided.
        
        :param list sample_ids:
            A list of sample IDs, e.g. ``[('A', 1), ('A', 2), ...]``.
        :param bool process:
            Process ``sample_ids`` by the ``sample_id_processor`` method.
        :param bool return_idx:
            Return the index array corresponding to the sort.
        """
        
        idx = self.argsort_by_sample_id(
            sample_ids = sample_ids,
            process = process,
        )
        
        self.sort_by_index(idx)
        
        if return_idx:
            
            return idx
    
    def sort_by(self, other, return_idx = False):
        """
        Sorts the sample attributes according to an other ``SampleSetAttrs``
        object. The other object must have the same sample IDs.
        
        :param SampleSetAttrs other:
            An other ``SampleSetAttrs`` with the same sample IDs.
        :param bool return_idx:
            Return the index array corresponding to the sort.
        """
        
        return self.sort_by_sample_id(
            other.sample_index_to_id,
            return_idx = return_idx,
        )


class SampleSorter(object):
    
    def __init__(
            self,
            sample_data = None,
            sample_ids = None,
            sample_id_processor = None,
            sample_axis = 0,
        ):
        """
        Keeps the order of samples synchronized between multiple objects.
        These objects represent sets of the same samples such as
        ``sample.SampleSet`` or ``feature.SampleData``.
        
        :param list,set sample_data:
            Other ``sample.SampleSet`` or ``feature.SampleData`` derived
            objects that should keep the same order of samples.
        :param list sample_ids:
            A list of sample IDs which determines the initial ordering. All
            objects in ``sample_data`` will be sorted according to this order.
            If not provided the ordering in the first object in
            ``sample_data`` will be used. You must provide at least
            ``sample_ids`` or ``sample_data``.
        :param int sample_axis:
            Which axis in the arrays corresponds to the samples.
            In ``sample.SampleSet`` objects this is axis 1 as axis 0
            corresponds to the features. In ``feature.SampleData`` derived
            objects this is axis 0.
        """
        
        self._sample_data = {}
        self._sample_axis = sample_axis
        
        if sample_data is None:
            
            sample_data = []
        
        if not isinstance(sample_data, (list, set)):
            
            sample_data = [sample_data]
        
        if not hasattr(self, 'attrs'):
            
            self._init_attrs(
                sample_data = sample_data,
                sample_ids  = sample_ids,
                sample_id_processor = sample_id_processor,
            )
        
        for s in sample_data:
            
            self.register(s)
    
    def _init_attrs(self, sample_data, sample_ids, sample_id_processor):
        """
        Initializes a ``SampleSetAttrs`` object to keep track of the order
        of samples by their IDs and ensure all objects will be set to the
        same ordering at initialization.
        """
        
        if sample_ids:
            
            self.attrs = SampleSetAttrs(
                sample_ids = sample_ids,
                sample_id_processor = sample_id_processor,
            )
            
        elif sample_data:
            
            first = sample_data[0]
            
            self.attrs = SampleSetAttrs(
                sample_ids = first.attrs.sample_id_to_index,
                sample_id_processor = first.attrs._sample_id_processor,
            )
            
        else:
            
            raise RuntimeError(
                'SampleSorter: `sample_data` or `sample_ids` '
                'must be provided.'
            )
    
    def sort_by(self, s):
        """
        Makes sure object ``s`` has the same ordering as all in this sorter.
        """
        
        idx = s.attrs.sort_by(self.attrs, return_idx = True)
        
        s._sort(idx)
    
    def _sort(self, idx):
        """
        Sorts only variables in this object by indices.
        """
        
        if hasattr(self, 'var'):
            
            numof_samples = len(self.attrs)
            
            for var in self.var:
                
                arr = getattr(self, var)
                
                if (
                    len(arr.shape) <= self._sample_axis or
                    arr.shape[self._sample_axis] != numof_samples
                ):
                    
                    continue
                
                setattr(
                    self,
                    var,
                    np.take(arr, idx, axis = self._sample_axis),
                )
    
    def register(self, s):
        """
        Registers a ``sample.SampleSet`` or ``feature.SampleData`` derived
        object ensuring it will keep the same order of samples.
        
        :param SampleSet,SampleData s:
            A ``sample.SampleSet`` or ``feature.SampleData`` derived
            object.
        """
        
        self._sort_by(s)
        
        self._sample_data[id(s)] = s
        
        if id(self) not in s._sample_data:
            
            s.register(self)
    
    def sort_by_sample_ids(self, sample_ids, process = False):
        """
        Sorts all connected objects by a list of sample IDs.
        
        :param list sample_ids:
            A list of sample ids in the desired order, e.g.
            ``[('A', 1), ('A', 2), ...]``.
        :param bool process:
            Use the ``sample_id_processor`` methods to convert the sample IDs
            provided.
        """
        
        idx = self.attrs.argsort_by_sample_id(
            sample_ids,
            process = process,
        )
        
        self.sort(idx)
    
    def sort(self, idx, _done = None):
        """
        Sorts all connected objects by indices.
        
        :param list idx:
            A list of indices.
        :param set _done:
            As the sorting propagates across objects this ``set`` keeps track
            which objects have been already sorted. Should be ``None`` when
            called by user.
        """
        
        _done = set() if _done is None else _done
        
        if id(self) in _done:
            
            return
        
        numof_samples = self.numof_samples
        
        if len(idx) != numof_samples:
            
            raise RuntimeError(
                'Invalid index length: %u while number of samples is %u.' % (
                    len(idx), numof_samples
                )
            )
        
        self.attrs.sort_by_index(idx)
        self._sort(idx)
        
        _done.add(id(self))
        
        for sd_id, sd in iteritems(self._sample_data):
            
            if sd_id not in _done:
                
                sd.order_samples(idx = idx, _done = _done)
