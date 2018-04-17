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
    
    def __init__(self, c = (2, 36), u = (0, 10), counts = {'H': -2}, **kwargs):
        
        metabolite.AbstractSubstituent.__init__(
            self, cores = ['O'], counts = counts, c = c, u = u, **kwargs
        )


class HydroxyFattyAcyl(metabolite.AbstractSubstituent):
    
    def __init__(
            self,
            c = (2, 24),
            u = (0, 6),
            counts = {'H': -2, 'O': 1},
            **kwargs
        ):
        
        metabolite.AbstractSubstituent.__init__(
            self, cores = ['O'],
            counts = counts,
            c = c,
            u = u,
            prefix = '2OH', **kwargs
        )


class FattyAlkoxy(metabolite.AbstractSubstituent):
    
    def __init__(self, c = (2, 36), u = (0, 10), counts = {}, **kwargs):
        
        metabolite.AbstractSubstituent.__init__(
            self, cores = [''], counts = counts, c = c, u = u, **kwargs
        )


class Sphingosine(metabolite.AbstractSubstituent):
    
    def __init__(
            self,
            c = (12, 24),
            u = (1, 6),
            counts = {},
            prefix = 'd',
            **kwargs
        ):
        """
        Note: this is a sphingosine backbone in a sphingolipid
        hence 2 hydrogens are missing and these should be replaced
        by the appropriate substituents.
        """
        
        if u[0] == 0:
            
            u = (1, u[1])
        
        metabolite.AbstractSubstituent.__init__(
            self, cores = ['O2N'], counts = counts, c = c, u = u,
            prefix = prefix, **kwargs
        )
    
    def get_prefix(self):
        
        return 'DH' if self.u == 0 and self.prefix == 'd' else self.prefix


class DihydroSphingosine(Sphingosine):
    
    def __init__(self, c = (12, 24), u = (0, 6), counts = {}, **kwargs):
        
        Sphingosine.__init__(
            self, c = c, u = u, counts = counts, prefix = 'DH', **kwargs
        )


class HydroxySphingosine(Sphingosine):
    
    def __init__(self, c = (12, 24), u = (0, 6), counts = {'O': 1}):
        
        Sphingosine.__init__(
            self, c = c, u = u, counts = counts, prefix = 't'
        )
