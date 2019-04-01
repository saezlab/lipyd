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


def ppm_tolerance(ppm, m):
    """
    Converts ppm value to plus/minus range.

    Parameters
    ----------
    int : ppm
        Tolerance in ppm.
    float : m
        The m/z value.

    Returns
    -------
    The m/z difference corresponding to the ppm value at this particular m/z.
    """
    
    return m * ppm * 1e-6


def findall(a, m, t = 20):
    """
    Finds all values within a given range of tolerance around a reference
    value in a one dimensional sorted numpy array (-slice) of floats.
    Returns list of array indices.

    Parameters
    ----------
    a : numpy
        Sorted one dimensional float array.
    m : float
        Value to lookup.
    t : float
        Range of tolerance (highest accepted difference) in ppm.
        (Default value = 20)

    Returns
    -------
    A list of array indices.
    """
    
    t_abs = ppm_tolerance(t, m)
    
    return _findall(a, m, t_abs)


def _findall(a, m, t):
    
    result = []
    
    # the upper closest index
    iu = a.searchsorted(m)
    
    if a.shape[0] > iu:
        
        u = 0
        
        while True:
            
            if iu + u < a.shape[0] and a[iu + u] - m <= t:
                
                result.append(iu + u)
                u += 1
                
            else:
                
                break
    
    if iu > 0:
        
        l = 1
        
        while True:
            
            if iu - l >= 0 and m - a[iu - l] <= t:
                
                result.append(iu - l)
                l += 1
                
            else:
                
                break
    
    return result


def find(a, m, t = 20):
    """Finds closest value based on a reference value in a one dimensional
    sorted numpy array of floats.
    Returns the array index closest value or `None` if nothing found
    in the range of tolerance.
    If the array contains contains more identical elements only the
    index of the first is returned.
    To lookup all values in a certain range see `lookup.findall()`.

    Parameters
    ----------
    a : numpy.array
        Sorted one dimensional float array (-slice).
    m : float
        Value to lookup.
    t : float
        Range of tolerance (highest accepted difference) in ppm.


    Returns
    -------
    Index of the closest value, None if no value found in within tolerance.
    """
    
    t_abs = ppm_tolerance(t, m)
    
    return _find(a, m, t_abs)


def _find(a, m, t):
    
    iu = a.searchsorted(m)
    
    dl = du = 9999.
    
    if iu < len(a):
        
        du = abs(a[iu] - m)
    
    if iu != 0:
        
        dl = abs(m - a[iu - 1])
    
    if dl < du:
        
        if dl < t:
            
            return iu - 1
    
    elif du <= t and iu < len(a):
        
        return iu


def match(observed, theoretical, tolerance = 20):
    """

    Parameters
    ----------
    observed :
        
    theoretical :
        
    tolerance :
         (Default value = 20)

    Returns
    -------
    Boolean.
    """
    
    tolerance_abs = ppm_tolerance(tolerance, theoretical)
    
    return _match(observed, theoretical, tolerance_abs)


def _match(observed, theoretical, tolerance):
    
    return abs(theoretical - observed) <= tolerance
