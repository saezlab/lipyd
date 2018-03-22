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


class FattyAcyl(metabolite.AbstractSubstituent):
    
    def __init__(self, c = (1, 25), u = (0, 7), counrs = {'H': -2}):
        
        metabolite.AbstractSubstituent.__init__(
            self, cores = ['O2'], counts = {'H': -2}, c = c, u = u
        )

class AbstractGlycerol(metabolite.AbstractMetabolite):
    
    def __init__(self, headgroup_weight = None,
                 positive = 0, negative = 0,
                 name = 'Lipid', num_fa = 2, type = 'GPL'):
        


class AbstractGPL(AbstractGlycerol):
    
    def __init__(self, ):
        
        AbstractGlycerol.__init__(self, c2 = , c3 = )
        self.bind('PO4')
