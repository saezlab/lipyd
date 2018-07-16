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

import lipyd.metabolite as metabolite


class FattyAcyl(metabolite.AbstractSubstituent):
    
    def __init__(self, c = (2, 36), u = (0, 10), counts = None, **kwargs):
        
        metabolite.AbstractSubstituent.__init__(
            self,
            cores = ['O'],
            counts = counts or {'H': -2},
            c = c,
            u = u,
            **kwargs
        )


class HydroxyFattyAcyl(metabolite.AbstractSubstituent):
    
    def __init__(
            self,
            c = (2, 24),
            u = (0, 6),
            counts = None,
            **kwargs
        ):
        
        metabolite.AbstractSubstituent.__init__(
            self,
            cores = ['O'],
            counts = counts or {'H': -2, 'O': 1},
            c = c,
            u = u,
            prefix = '2OH',
            **kwargs
        )


class FattyAlkoxy(metabolite.AbstractSubstituent):
    
    def __init__(self, c = (2, 36), u = (0, 10), counts = None, **kwargs):
        
        metabolite.AbstractSubstituent.__init__(
            self,
            cores = [''],
            counts = counts or {},
            c = c,
            u = u,
            prefix = 'O-',
            **kwargs
        )


class Sphingosine(metabolite.AbstractSubstituent):
    
    def __init__(
            self,
            c = (12, 24),
            u = (1, 6),
            counts = None,
            prefix = 'd',
            keto = False,
            **kwargs
        ):
        """
        Note: this is a sphingosine backbone in a sphingolipid
        hence 2 hydrogens are missing and these should be replaced
        by the appropriate substituents.
        """
        
        counts = counts or {}
        
        if keto:
            counts['H'] = counts['H'] - 2 if 'H' in counts else -2
        
        metabolite.AbstractSubstituent.__init__(
            self, cores = ['O2N'], counts = counts, c = c, u = u,
            prefix = prefix, **kwargs
        )
    
    def get_prefix(self):
        
        return 'DH' if self.u == 0 and self.prefix == 'd' else self.prefix


class DihydroSphingosine(Sphingosine):
    
    def __init__(self, c = (12, 24), u = (0, 6), counts = None, **kwargs):
        
        Sphingosine.__init__(
            self, c = c, u = u, counts = counts or {}, prefix = 'DH', **kwargs
        )


class HydroxySphingosine(Sphingosine):
    
    def __init__(self, c = (12, 24), u = (0, 6), counts = None):
        
        Sphingosine.__init__(
            self, c = c, u = u, counts = counts or {'O': 1}, prefix = 't'
        )
