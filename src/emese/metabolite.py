#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `emese` python module
#
#  Copyright (c) 2015-2017 - EMBL
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
from past.builtins import xrange, range

import itertools

import emese.mass as mass

class AbstractMetabolite(mass.Formula):
    
    def __init__(self, core = None,
                 charge = 0,
                 isotope = 0,
                 name = 'Unknown',
                 syn = {},
                 subs = [],
                 **kwargs):
        
        mass.MolWeight.__init__(self,
            core if type(core) is not float else None,
            charge,
            isotope,
            **kwargs)
        
        if not self.has_weight() and type(core) is float:
            
            self.weight = core
            
        else:
            
            raise ValueError('Please provide either formula or '
                             'atom counts or weight.')
        
        self.name = name
        self.syn = syn
        self.subs = subs


class MetaboliteClass(AbstractMetabolite):
    
    def __init__(self):
        
        
        

class SubstituentSeries(object):
    
    def __init__(self, default = 'H', variations = set([]),
                 charge = 0, isotope = 0):
        
        self.charge  = self.get_range(charge)
        self.isotope = self.get_range(isotope)
        
        
    
    @staticmethod
    def get_range(i):
        
        return i if type(i) is tuple else (i, i)


class Metabolite(MetaboliteClass):
    
    def __init__(self):
        
        
