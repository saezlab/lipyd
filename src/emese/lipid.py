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

import emese.mass as mass

class AbstractMetabolite(mass.MolWeight):
    
    def __init__(self, formula = None, weight = None,
                 charge = 0, isotope = 0,
                 name = 'Unknown', synonyms = {}, **kwargs):
        
        mass.MolWeight.__init__(self, formula, charge, isotope, **kwargs)
        
        if not self.has_weight() and weight:
            
            self.weight = weight
        
        self.name = name
        self.synonyms = synonyms


class MetaboliteClass(AbstractMetabolite):
    
    def __init__(self):
        
        


class Metabolite(MetaboliteClass):
    
    def __init__(self):
        
        


class Lipid(Metabolite):
    
    def __init__(self, headgroup_weight = None,
                 positive = 0, negative = 0,
                 name = 'Lipid', num_fa = 2, type = 'GPL'):
        
        
