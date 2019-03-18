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

"""
Module for arithmetics with m/z values and adduct ions.
"""

from __future__ import print_function
from future.utils import iteritems
from past.builtins import xrange, range, reduce

import imp

import lipyd.mass as mass

class Mz():
    """Represents one m/z value.
    Provides methods for conversion to exact mass or various adducts.
    Other classes representing molecular entities and masses inherit
    these methods from here.

    Parameters
    ----------

    Returns
    -------

    """
    
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
        """

        Parameters
        ----------
        m :
            

        Returns
        -------

        """
        return (self.mz * self.z + float(m)) / abs(self.z)
    
    def weight(self):
        """ """
        return self.mz * self.z
    
    def remove_h(self):
        """ """
        return self.adduct(-mass.proton)
    
    def remove_ac(self):
        """ """
        m = mass.MassBase('H3C2O2')
        return self.adduct(-m - mass.electron)
    
    def remove_fo(self):
        """ """
        m = mass.MassBase('HCO2')
        return self.adduct(-m - mass.electron)
    
    def remove_nh4(self):
        """ """
        m = mass.MassBase('NH4')
        return self.adduct(-m + mass.electron)
    
    def remove_oh(self):
        """ """
        m = mass.MassBase('OH')
        return self.adduct(-m - mass.electron)
    
    def add_h(self):
        """ """
        return self.adduct(mass.proton)
    
    def add_2h(self):
        """ """
        return self.adduct(2 * mass.proton)
    
    def add_3h(self):
        """ """
        return self.adduct(3 * mass.proton)
    
    def add_oh(self):
        """ """
        m = mass.MassBase('OH')
        return self.adduct(m + mass.electron)
    
    def add_fo(self):
        """ """
        m = mass.MassBase('HCO2')
        return self.adduct(m + mass.electron)
    
    def add_fo_nafo(self):
        """ """
        m = mass.MassBase('HCO2NaHCO2')
        return self.adduct(m + mass.electron)
    
    def remove_fo_nafo(self):
        """ """
        m = mass.MassBase('HCO2NaHCO2')
        return self.adduct(-m - mass.electron)
    
    def add_fo_2nafo(self):
        
        m = mass.MassBase('HCO2Na2H2C2O4')
        return self.adduct(m + mass.electron)
    
    def remove_fo_2nafo(self):
        
        m = mass.MassBase('HCO2Na2H2C2O4')
        return self.adduct(-m - mass.electron)
    
    def add_fo_3nafo(self):
        
        m = mass.MassBase('HCO2Na3H3C3O6')
        return self.adduct(m + mass.electron)
    
    def remove_fo_3nafo(self):
        
        m = mass.MassBase('HCO2Na3H3C3O6')
        return self.adduct(-m - mass.electron)
    
    def add_fo_4nafo(self):
        
        m = mass.MassBase('HCO2Na4H4C4O8')
        return self.adduct(m + mass.electron)
    
    def remove_fo_4nafo(self):
        
        m = mass.MassBase('HCO2Na4H4C4O8')
        return self.adduct(-m - mass.electron)
    
    def add_fo_5nafo(self):
        
        m = mass.MassBase('HCO2Na5H5C5O10')
        return self.adduct(m + mass.electron)
    
    def remove_fo_5nafo(self):
        
        m = mass.MassBase('HCO2Na5H5C5O10')
        return self.adduct(-m - mass.electron)
    
    def remove_h_add_nafo(self):
        
        m = mass.MassBase('NaCO2')
        return self.adduct(m + mass.electron)
    
    def add_h_remove_nafo(self):
        
        m = mass.MassBase('NaCO2')
        return self.adduct(-m - mass.electron)
    
    def remove_h_add_2nafo(self):
        
        m = mass.MassBase('Na2HC2O4')
        return self.adduct(m + mass.electron)
    
    def add_h_remove_2nafo(self):
        
        m = mass.MassBase('Na2HC2O4')
        return self.adduct(-m - mass.electron)
    
    def remove_h_add_3nafo(self):
        
        m = mass.MassBase('Na3H2C3O6')
        return self.adduct(m + mass.electron)
    
    def add_h_remove_3nafo(self):
        
        m = mass.MassBase('Na3H2C3O6')
        return self.adduct(-m - mass.electron)
    
    def remove_h_add_na_3nafo(self):
        
        m = mass.MassBase('Na4HC3O6')
        return self.adduct(m + mass.electron)
    
    def add_h_remove_na_3nafo(self):
        
        m = mass.MassBase('Na4HC3O6')
        return self.adduct(-m - mass.electron)
    
    def remove_h_add_na_4nafo(self):
        
        m = mass.MassBase('Na5H2C4O8')
        return self.adduct(m + mass.electron)
    
    def add_h_remove_na_4nafo(self):
        
        m = mass.MassBase('Na5H2C4O8')
        return self.adduct(-m - mass.electron)
    
    def remove_h_add_na_5nafo(self):
        
        m = mass.MassBase('Na6H3C5O10')
        return self.adduct(m + mass.electron)
    
    def add_h_remove_na_5nafo(self):
        
        m = mass.MassBase('Na6H3C5O10')
        return self.adduct(-m - mass.electron)
    
    def remove_h_add_na_6nafo(self):
        
        m = mass.MassBase('Na7H4C6O12')
        return self.adduct(m + mass.electron)
    
    def add_h_remove_na_6nafo(self):
        
        m = mass.MassBase('Na7H4C6O12')
        return self.adduct(-m - mass.electron)
    
    def remove_h_add_na_7nafo(self):
        
        m = mass.MassBase('Na8H5C7O14')
        return self.adduct(m + mass.electron)
    
    def add_h_remove_na_7nafo(self):
        
        m = mass.MassBase('Na8H5C7O14')
        return self.adduct(-m - mass.electron)
    
    def add_ac(self):
        """ """
        m = mass.MassBase('H3C2O2')
        return self.adduct(m + mass.electron)
    
    def add_nh4(self):
        """ """
        m = mass.MassBase('NH4')
        return self.adduct(m - mass.electron)
    
    def add_na(self):
        """ """
        m = mass.MassBase('Na')
        return self.adduct(m - mass.electron)
    
    def remove_na(self):
        """ """
        m = mass.MassBase('Na')
        return self.adduct(-m + mass.electron)
    
    def reload(self):
        """ """
        modname = self.__class__.__module__
        mod = __import__(modname, fromlist=[modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)
