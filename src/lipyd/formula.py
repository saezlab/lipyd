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

from future.utils import iteritems

from collections import defaultdict
from argparse import Namespace
import copy

import lipyd.mass as mass
import lipyd.mz as mz


def formula2atoms(formula):
    """
    Converts chemical formula string to dict of atom counts.
    
    Args
    ----
    :param str formula:
        Chemical formula
    """
    
    atoms = defaultdict(lambda: 0)
    
    for elem, cnt in mass.reform.findall(formula):
        
        atoms[elem] += int(cnt or '1')
    
    return atoms


class Formula(mass.MassBase, mz.Mz):
    
    def __init__(
            self,
            formula = None,
            charge = 0,
            isotope = 0,
            z = 1,
            sign = None,
            tolerance = .01,
            attrs = None,
            **kwargs
        ):
        
        attrs = attrs or {}
        attrs.update(kwargs)
        self._set_attrs(attrs)
        
        if isinstance(formula, Formula):
            
            _attrs = copy.deepcopy(formula.attrs.__dict__)
            charge = formula.charge
            isotope = formula.isotope
            formula = formula.formula
            self.attrs = _attrs.update(self.attrs.__dict__)
        
        mass.MassBase.__init__(self, formula, charge, isotope, **kwargs)
        
        if self.formula == '':
            
            self.mass = 0.0
            self.mass_calculated = True
            
        elif self.mass == 0.0:
            
            self.formula = ''
            self.mass_calculated = True
        
        self.reset_atoms()
        self.add(self.formula if self.formula else '')
        
        mz.Mz.__init__(
            self,
            mz = self.mass / z,
            z = z,
            sign = sign,
            tolerance = tolerance,
        )
    
    def __add__(self, other):
        
        if (
            not self.has_formula() or
            type(other) is float or
            (
                hasattr(other, 'has_formula') and
                not other.has_formula()
            )
        ):
            
            new_mass = mass.MassBase.__add__(self, other)
            new_charge = self.charge + (
                other.charge
                if hasattr(other, 'charge')
                else 0
            )
            new_isotope = self.isotope + (
                other.isotope
                if hasattr(other, 'isotope')
                else 0
            )
            new = Formula(
                new_mass,
                charge = new_charge,
                isotope = new_isotope,
                attrs = self.attrs.__dict__
            )
            
        else:
            
            new = Formula('%s%s' % (
                    self.formula,
                    other.formula
                        if hasattr(other, 'formula')
                        else other
                ),
                charge = self.charge + (
                    other.charge
                    if hasattr(other, 'charge')
                    else 0
                ),
                isotope = self.isotope + (
                    other.isotope
                    if hasattr(other, 'isotope')
                    else 0
                ),
                attrs = self.attrs.__dict__
            )
            
            for a in ('c', 'u'):
                if a in self.attrs and a in other.attrs:
                    setattr(
                        new.attrs,
                        a,
                        getattr(self.attrs, a) + getattr(other.attrs, a)
                    )
        
        if new.mass == 0.0 or new.formula == '':
            
            new.formula = ''
            new.mass_calculated = True
        
        new.update_mz(
            z = self.z,
            sign = self.sign,
            tolerance = self.tol
        )
        
        new._set_attrs
        
        return new
    
    def __iadd__(self, other):
        
        if type(other) is float:
            
            self.mass = self.mass + other
            self.formula = ''
            self.reset_atoms()
            self.mass_calculated = False
            
        elif not self.has_formula():
            
            self.mass += other.mass
            
        else:
            
            self.add(other.formula if hasattr(other, 'formula') else other)
        
        self.charge += (other.charge if hasattr(other, 'charge') else 0)
        self.isotope += (other.isotope if hasattr(other, 'isotope') else 0)
        
        self.calc_mass()
        self.update_mz()
        
        return self
    
    def __sub__(self, other):
        
        new = copy.copy(self)
        new.__isub__(other)
        return new
    
    def __isub__(self, other):
        
        if type(other) is float:
            
            self.mass = self.mass - other
            self.formula = ''
            self.reset_atoms()
            self.mass_calculated = False
            
        elif not self.has_formula():
            
            self.mass -= other.mass
            
        else:
            
            self.sub(other.formula if hasattr(other, 'formula') else other)
        
        self.charge -= (other.charge if hasattr(other, 'charge') else 0)
        self.isotope -= (other.isotope if hasattr(other, 'isotope') else 0)
        
        return self
    
    def __imul__(self, other):
        
        if not type(other) is int:
            
            return self
        
        for elem, cnt in iteritems(self.atoms):
            
            self.counts[elem] = cnt * other
        
        self.calc_mass()
        self.formula_from_dict(self.atoms)
        self.isotope = self.isotope * other
    
    def __mul__(self, other):
        
        if not type(other) is int:
            
            return copy.deepcopy(self)
        
        new_atoms = defaultdict(int)
        
        for elem, cnt in iteritems(self.atoms):
            
            new_atoms[elem] = cnt * other
        
        return Formula(
            **new_atoms,
            isotope = self.isotope * other,
            charge = self.charge
        )
    
    def __iter__(self, **kwargs):
        
        yield self
    
    def reset_atoms(self):
        
        self.atoms = defaultdict(int)
    
    def as_mass(self):
        
        return MassBase(self.formula, self.charge, self.isotope)
    
    def add(self, formula):
        
        for elem, cnt in mass.reform.findall(formula):
            self.atoms[elem] += int(cnt or '1')
        
        self.update()
    
    def sub(self, formula):
        
        for elem, cnt in mass.reform.findall(formula):
            self.atoms[elem] -= int(cnt or '1')
            
            if self.atoms[elem] < 0:
                
                raise ValueError('Can not remove %s from %s: '
                    'too few %s atoms!' % (formula, self.formula, elem))
        
        self.update()
    
    def update(self):
        
        if len(self.atoms):
            
            self.formula = ''.join('%s%u' % (elem, self.atoms[elem])
                                    for elem in sorted(self.atoms.keys()))
            self.calc_mass()
    
    def bind(self, other, loss = 'H2O'):
        
        return self + other - loss
    
    def split(self, product1, add = 'H2O'):
        
        product1 = Formula(product1)
        
        return product1, self - product1 + add
    
    def update_mz(self, mz = None, z = 1, sign = None, tolerance = .01,
                  overwrite = False):
        
        if not hasattr(self, 'z') or overwrite:
            self.z = z
        
        if not hasattr(self, 'sign') or overwrite:
            self.sign = sign
        
        if not hasattr(self, 'tolerance') or overwrite:
            self.tol = tolerance
        
        self.mz = mz or self.mass / self.z
    
    def getname(self):
        
        return self.name if hasattr(self, 'name') else self.__str__()
    
    def __str__(self):
        
        return str(self.formula)
    
    def _set_attrs(self, attrs):
        
        attrs = attrs or {}
        self.attrs = (
            attrs
            if hasattr(attrs, '__dict__') else
            Namespace(**attrs)
        )


class Mass(Formula):
    
    def __init__(self, formula_mass = None, charge = 0, isotope = 0, **kwargs):
        
        if ((formula_mass is None and not kwargs) and
            type(formula_mass) is float):
            
            # unknown formula, initializing an empty Formula:
            Formula.__init__(self, '', charge = charge, isotope = isotope)
            self.mass = formula_mass
            
        else:
            
            Formula.__init__(self, formula_mass,
                             charge = charge,
                             isotope = isotope,
                             **kwargs)
    
    def bind(self, other, loss = 'H2O'):
        
        if self.has_formula() and (type(other) is str):
            
            pass
