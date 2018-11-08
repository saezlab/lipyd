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

import re
import numpy as np

IONMODE_POS = 'pos'
IONMODE_NEG = 'neg'

refloat = re.compile(r'([-]?[0-9]*[\.]?[0-9]+[eE]?[-\+]?[0-9]*)')
reint   = re.compile(r'([-]?[0-9]+[\.]?[0-9]*)')

try:
    basestring
except NameError:
    basestring = str

def guess_ionmode(*args):
    
    for a in args:
        
        if hasattr(a, 'lower'):
            
            a = a.lower()
            
            if IONMODE_POS in a:
                
                return IONMODE_POS
                
            elif IONMODE_NEG in a:
                
                return IONMODE_NEG


def to_float(num):
    """
    Extracts ``float`` from string, or returns ``numpy.nan``.
    """
    
    if isinstance(num, float):
        
        return num
    
    if isinstance(num, int):
        
        return float(num)
    
    if isinstance(num, basestring):
        
        num = num.strip()
        match = refloat.match(num)
        
        if match:
            
            return float(match.groups()[0])
            
        else:
            
            if num.lower() == 'inf':
                
                return np.inf
            
            if num.lower() == '-inf':
                
                return -np.inf
    
    return np.nan


def to_int(num):
    """
    Extracts ``int`` from string.
    """
    
    if isinstance(num, int):
        
        return num
    
    match = reint.match(num.strip())
    
    if match:
        
        return int(match.groups(0)[0])
        
    else:
        
        raise ValueError('Integer expected: %g' % num)
