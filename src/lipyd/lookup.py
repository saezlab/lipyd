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


def ppm_tolerance(ppm, m):
    """
    Converts ppm value to plus/minus range.
    
    :param int ppm:
        Tolerance in ppm.
    :param float m:
        The m/z value.
    """
    
    return m * ppm * 1e-6


def findall(a, m, t = 20):
    """
    Finds all values within a given range of tolerance around a reference
    value in a one dimensional sorted numpy array (-slice) of floats.
    Returns list of array indices.
    
    :param numpy.array a: Sorted one dimensional float array.
    :param float m: Value to lookup.
    :param float t: Range of tolerance (highest accepted difference) in ppm.
    """
    
    result = []
    
    t = ppm_tolerance(t, m)
    
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
    """
    Finds closest value based on a reference value in a one dimensional
    sorted numpy array of floats.
    Returns the array index closest value or `None` if nothing found
    in the range of tolerance.
    If the array contains contains more identical elements only the
    index of the first is returned.
    To lookup all values in a certain range see `lookup.findall()`.
    
    :param numpy.array a: Sorted one dimensional float array (-slice).
    :param float m: Value to lookup.
    :param float t: Range of tolerance (highest accepted difference) in ppm.
    """
    
    t = ppm_tolerance(t, m)
    
    iu = a.searchsorted(m)
    
    dl = du = None
    
    du = abs(a[iu] - m)
    
    if iu != 0:
        
        dl = abs(m - a[iu - 1])
    
    if dl is not None and dl < du:
        
        if dl < t:
            
            return iu - 1
    
    elif du <= t:
        
        return iu


def match(observed, theoretical, tolerance = 20):
    
    tolerance = ppm_tolerance(tolerance, theoretical)
    
    return abs(theoretical - observed) <= tolerance
