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
    
    def __init__(self):
        
        
