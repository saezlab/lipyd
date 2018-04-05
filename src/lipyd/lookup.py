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

def findall(a, m, t = .01):
    """
    Finds all values within a given range of tolerance around a reference
    value in a one dimensional sorted numpy array of floats.
    Returns list of floats.
    
    :param numpy.array a: Sorted one dimensional float array.
    :param float m: Value to lookup.
    :param float t: Range of tolerance (highest accepted difference).
    """
    
    result = []
    
    # the upper closest index
    iu = a.searchsorted(m)
    
    if a.shape[0] > iu:
        
        u = 0
        
        while True:
            
            if a[iu + u] - m <= t:
                
                result.append(a[iu + u])
                u += 1
                
            else:
                
                break
    
    if iu > 0:
        
        l = 1
        
        while True:
            
            if iu - l >= 0 and m - a[iu - l] <= t:
                
                result.append(a[iu - l])
                l += 1
                
            else:
                
                break
    
    return result


def find(a, m, t = .01):
    """
    Finds closest value based on a reference value in a one dimensional
    sorted numpy array of floats.
    Returns closest value or `None` if nothing found in the range of
    tolerance.
    To lookup all values in a certain range see `lookup.findall()`.
    
    :param numpy.array a: Sorted one dimensional float array.
    :param float m: Value to lookup.
    :param float t: Range of tolerance (highest accepted difference).
    """
    
    iu = a.searchsorted(m)
    
    dl = du = None
    
    du = abs(a[iu] - m)
    
    if iu != 0:
        
        dl = abs(m - a[iu - 1])
    
    if dl is not None and dl < du:
        
        if dl < t:
            
            return a[iu - 1]
    
    elif du <= t:
        
        return a[iu]
