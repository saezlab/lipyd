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

from past.builtins import xrange, range

import itertools

import numpy as np

import lipyd.lookup as lookup


def align(*arrays, tolerance = 5):
    
    narrays = len(arrays)
    
    array  = np.concatenate(arrays)
    # indices within arrays
    idx    = np.concatenate(tuple(
        np.arange(len(a))
        for a in arrays
    ))
    # which array the value belongs to
    iarray = np.concatenate(tuple(
        np.full(a.shape, i)
        for i, a in enumerate(arrays)
    ))
    # keep track if elements have been used already
    used = np.full(iarray.shape, False)
    
    # sort all arrays by values
    isort  = array.argsort()
    array  = array[isort]
    idx    = idx[isort]
    iarray = iarray[isort]
    
    # here we collect the lines of the result array
    result = []
    
    # these are the array memberships
    for i, val in enumerate(array):
        
        partners = lookup.findall(array, val)
        
        result.append((i, partners))
    
    return array, result
