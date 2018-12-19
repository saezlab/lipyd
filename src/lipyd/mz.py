#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `lipyd` python module
#
#  Copyright (c) 2015-2018 - EMBL
#
#  File author(s): Dénes Türei (turei.denes@gmail.com)
#
#  Distributed under the GNU GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://www.ebi.ac.uk/~denes
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
