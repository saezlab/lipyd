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

from __future__ import print_function
from future.utils import iteritems
from past.builtins import xrange, range, reduce

import imp

import emese.mass as mass

class Mz():
    
    def __init__(self, mz, z = 1, sign = None, tolerance = 0.01):
        self.mz = mz
        self.z = z
        self.sign = sign
        self.tol = tolerance
    
    def __eq__(self, other):
        return self.z == other.z and \
            self.mz > other.mz - self.tol and \
            self.mz < other.mz + self.tol
    
    def __str__(self):
        return 'm/z = %f' % self.mz
    
    def adduct(self, m):
        return (self.mz * self.z + float(m)) / abs(self.z)
    
    def weight(self):
        return self.mz * self.z
    
    def remove_h(self):
        return self.adduct(-mass.proton)
    
    def remove_ac(self):
        m = mass.Mass('H3C2O2')
        return self.adduct(-m - mass.electron)
    
    def remove_fo(self):
        m = mass.Mass('HCO2')
        return self.adduct(-m - mass.electron)
    
    def remove_nh4(self):
        m = mass.Mass('NH4')
        return self.adduct(-m + mass.electron)
    
    def remove_oh(self):
        m = mass.Mass('OH')
        return self.adduct(-m - mass.electron)
    
    def add_h(self):
        return self.adduct(mass.proton)
    
    def add_2h(self):
        return self.adduct(2 * mass.proton)
    
    def add_3h(self):
        return self.adduct(3 * mass.proton)
    
    def add_oh(self):
        m = mass.Mass('OH')
        return self.adduct(m + mass.electron)
    
    def add_fo(self):
        m = mass.Mass('HCO2')
        return self.adduct(m + mass.electron)
    
    def add_ac(self):
        m = mass.Mass('H3C2O2')
        return self.adduct(m + mass.electron)
    
    def add_nh4(self):
        m = mass.Mass('NH4')
        return self.adduct(m - mass.electron)
    
    def add_na(self):
        m = mass.Mass('Na')
        return self.adduct(m - mass.electron)
    
    def remove_na(self):
        m = mass.Mass('Na')
        return self.adduct(-m + mass.electron)
    
    def reload(self):
        modname = self.__class__.__module__
        mod = __import__(modname, fromlist=[modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)
