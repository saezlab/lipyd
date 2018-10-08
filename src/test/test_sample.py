#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `lipyd` python module
#
#  Copyright (c) 2015-2018 - EMBL
#
#  File author(s): Dénes Türei (turei.denes@gmail.com)
#
#  This code is not for public use.
#  Please do not redistribute.
#  For permission please contact me.
#
#  Website: http://www.ebi.ac.uk/~denes
#

import pytest

import warnings
import numpy as np

import lipyd.sample as sample


# here we test the sorting methods of the sample module with
# some mock data of length 7

mz = np.array([
    323.23,
    234.21,
    231.19,
    564.43,
    343.32,
    787.56,
    672.47,
])

intensitites = np.array([
    3000.0,
    12000.0,
    720000.0,
    30000.0,
    2430000.0,
    2900.0,
    4242000.0,
])


class TestSample(object):
    
    def test_feature_idx(self):
        
        idx = sample.FeatureIdx(10)
        
        assert np.all(idx._current == idx._original)
        
        rand10 = np.empty((10,))
        sort10 = rand10.argsort()
        idx.sort_all(sort10)
        
        assert np.all(idx._current == sort10)
        # the smallest number is at index zero currently
        # as arrays are sorted by the data
        assert idx._original[rand10.argmin()] == 0
        
        idx.sort_all()
        
        assert np.all(idx._current == idx._original)
    
    def test_feature_base(self):
        
        a0 = np.random.random(10)
        a1 = np.random.random(10)
        a2 = np.random.random(10)
        b1 = np.random.random(10)
        
        a0min = a0.argmin()
        a2a0min = a2[a0min]
        b1a0min = b1[a0min]
        
        f0 = sample.FeatureBase(a = a0)
        f1 = sample.FeatureBase(a = a1, b = b1, sorter = f0.sorter)
        f2 = sample.FeatureBase(a = a2, sorter = f0.sorter)
        
        f0.sort_all('a')
        
        assert f2.a[0] == a2a0min
        assert f1.b[0] == b1a0min
    
    def test_feature_base_3d(self):
        
        a0 = np.random.random((10, 7, 5))
        a1 = np.random.random((10, 7, 5))
        a2 = np.random.random((10, 7, 5))
        b1 = np.random.random((10, 7, 5))
        
        a0min = a0[:,0,0].argmin()
        a2a0min = a2[a0min,0,0]
        b1a0min = b1[a0min,0,0]
        
        f0 = sample.FeatureBase(a = a0)
        f1 = sample.FeatureBase(a = a1, b = b1, sorter = f0.sorter)
        f2 = sample.FeatureBase(a = a2, sorter = f0.sorter)
        
        f0.sort_all('a', indices = (0, 0))
        
        assert f2.a[0,0,0] == a2a0min
        assert f1.b[0,0,0] == b1a0min
        
        with pytest.warns(UserWarning):
            
            f2.sort_all('a')
