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

import lipyd.metabolite as metabolite
import lipyd.lipproc as lipproc


class FattyAcyl(metabolite.AbstractSubstituent):
    """ """
    
    def __init__(self, c = None, u = None, counts = None, **kwargs):
        
        c = c or (2, 36)
        u = u if u is not None else (0, 10)
        
        chain_attr = lipproc.ChainAttr()
        
        metabolite.AbstractSubstituent.__init__(
            self,
            cores = ['O'],
            counts = counts or {'H': -2},
            c = c,
            u = u,
            chain_attr = chain_attr,
            chain_type = 'FA',
            **kwargs,
        )


class HydroxyFattyAcyl(metabolite.AbstractSubstituent):
    """ """
    
    def __init__(
            self,
            c = None,
            u = None,
            counts = None,
            **kwargs
        ):
        
        c = c or (2, 24)
        u = u if u is not None else (0, 6)
        
        chain_attr = lipproc.ChainAttr(oh = ('OH',))
        
        metabolite.AbstractSubstituent.__init__(
            self,
            cores = ['O'],
            counts = counts or {'H': -2, 'O': 1},
            c = c,
            u = u,
            prefix = 'OH',
            chain_attr = chain_attr,
            chain_type = 'FAOH',
            **kwargs
        )


class FattyAlkoxy(metabolite.AbstractSubstituent):
    """ """
    
    def __init__(self, c = (2, 36), u = (0, 10), counts = None, **kwargs):
        
        chain_attr = lipproc.ChainAttr(ether = True)
        
        metabolite.AbstractSubstituent.__init__(
            self,
            cores = [''],
            counts = counts or {},
            c = c,
            u = u,
            prefix = 'O-',
            chain_attr = chain_attr,
            chain_type = 'FAL',
            **kwargs
        )


class Sphingosine(metabolite.AbstractSubstituent):
    """ """
    
    def __init__(
            self,
            c = None,
            u = None,
            counts = None,
            lcb_type = 'd',
            **kwargs
        ):
        """
        Note: this is a sphingosine backbone in a sphingolipid
        hence 2 hydrogens are missing and these should be replaced
        by the appropriate substituents.
        """
        
        c = c or (12, 24)
        u = u if u is not None else (1, 6)
        counts = counts or {}
        
        if lcb_type == 'k':
            counts['H'] = counts['H'] - 2 if 'H' in counts else -2
        if lcb_type == 't':
            counts['O'] = counts['O'] + 1 if 'O' in counts else 1
        
        chain_attr = lipproc.ChainAttr(sph = lcb_type)
        
        metabolite.AbstractSubstituent.__init__(
            self,
            cores = ['O2N'],
            counts = counts,
            c = c,
            u = u,
            chain_attr = chain_attr,
            chain_type = 'Sph',
            prefix = lcb_type,
            **kwargs
        )
    
    def get_prefix(self):
        """ """
        
        return 'DH' if self.u == 0 and self.prefix == 'd' else self.prefix


class DihydroSphingosine(Sphingosine):
    """ """
    
    def __init__(self, c = None, u = None, counts = None, **kwargs):
        
        c = c or (12, 24)
        u = u if u is not None else (0, 6)
        
        Sphingosine.__init__(
            self,
            c = c,
            u = u,
            counts = counts or {},
            prefix = 'DH',
            **kwargs
        )


class HydroxySphingosine(Sphingosine):
    """ """
    
    def __init__(self, c = None, u = None, counts = None, **kwargs):
        
        c = c or (12, 24)
        u = u if u is not None else (0, 6)
        
        Sphingosine.__init__(
            self,
            c = c,
            u = u,
            counts = counts or {'O': 1},
            prefix = 't',
            **kwargs
        )


class Sterol(metabolite.AbstractSubstituent):
    
    
    def __init__(self, c = None, u = None, counts = None, **kwargs):
        
        c = c or 0
        u = u or 0
        
        metabolite.AbstractSubstituent.__init__(
            self,
            cores = ['O'],
            counts = counts or {'H': -2},
            c = c,
            u = u,
            chain_attr = chain_attr,
            chain_type = 'St',
            **kwargs,
        )
