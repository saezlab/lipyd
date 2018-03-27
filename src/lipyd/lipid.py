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

from future.utils import iteritems

import lipyd.metabolite as metabolite
import lipyd.substituent as substituent
import lipyd.formula as formula


class AbstractGlycerol(metabolite.AbstractMetabolite):
    
    def __init__(
            self,
            sn1 = 'OH',
            sn2 = 'OH',
            sn3 = 'OH',
            pcharge = 0,
            ncharge = 0,
            name = 'Glycerol',
            typ  = 'G',
            **kwargs
        ):
        
        self.pcharge = pcharge
        self.ncharge = ncharge
        self.netcharge = self.pcharge - self.ncharge
        
        metabolite.AbstractMetabolite.__init__(
            self,
            core = 'C3H5',
            subs = [
                self.get_substituent(sn1),
                self.get_substituent(sn2),
                self.get_substituent(sn3)
            ],
            name = name,
            charge = self.netcharge,
            **kwargs
        )


class AbstractGPL(AbstractGlycerol):
    
    def __init__(
            self,
            headgroup = 'H',
            lyso = False,
            ether = False,
            fa_args = {},
            name = 'GPL',
            typ  = 'GPL',
            **kwargs
        ):
        
        AbstractGlycerol.__init__(
            self,
            sn1  = substituent.FattyAcyl(**fa_args),
            sn2  = 'OH' if lyso else substituent.FattyAcyl(**fa_args),
            sn3  = 'PO4H%s' % headgroup,
            name = name,
            **kwargs
        )


class Phosphatidylethanolamine(AbstractGPL):
    
    def __init__(self, **kwargs):
        
        AbstractGPL.__init__(
            self,
            headgroup = 'C2H4NH2',
            name = 'PE',
            typ  = 'PE',
            **kwargs
        )


class AbstractSphingolipid(metabolite.AbstractMetabolite):
    
    def __init__(
            self,
            shp = None,
            o = 'H',
            n = 'H',
            pcharge = 0,
            ncharge = 0,
            name = 'Sphingolipid',
            typ  = 'SL',
            **kwargs
        ):
        
        self.pcharge = pcharge
        self.ncharge = ncharge
        self.netcharge = self.pcharge - self.ncharge
        
        metabolite.AbstractMetabolite.__init__(
            self,
            core = 'H',
            subs = [
                self.get_substituent(sn1),
                self.get_substituent(sn2),
                self.get_substituent(sn3)
            ],
            name = name,
            charge = self.netcharge,
            **kwargs
        )
