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

from past.builtins import xrange, range
from future.utils import iteritems

import bs4
import re
import warnings
import imp
import sys
import copy
import functools
import operator
import itertools
import collections

try:
    import pyopenms as oms
except:
    pass

import lipyd._curl as _curl


#: Mass of a proton
proton = 1.00727646677
#: Mass of an electron
electron = 0.00054857990924
#: Mass of a neutron
neutron = 1.00866491588

#: Mass of a proton
p = proton
#: Mass of an electron
e = electron
#: Mass of a neutron
n = neutron

_re_nondigit = re.compile(r'[^\d.]+')
_re_form  = re.compile(r'([A-Za-z][a-z]*)([0-9]*)')
replmi  = re.compile(r'([-+])')
refloat = re.compile(r'[0-9\.]+')


def formula_to_atoms(formula):
    """
    Converts chemical formula string to dict of atom counts.
    
    Parameters
    ----------
    formula : str
        Chemical formula, e.g. ``CH3COOH``.

    Returns
    -------
    ``dict`` with elements as keys and counts as values.
    """
    
    atoms = collections.defaultdict(int)
    
    for elem, cnt in _re_form.findall(formula):
        
        atoms[elem] += int(cnt or '1')
    
    return atoms


class MassDatabase(object):
    """
    Downloads, processes and serves data for atomic and isotopic masses
    and weights.
    """
    
    #: URL for atomic masses
    url_masses = 'http://www.ciaaw.org/atomic-masses.htm'
    #: URL for atomic weights
    url_weights = 'http://www.ciaaw.org/atomic-weights.htm'
    #: URL for isotopic abundances
    url_abundances = 'http://www.ciaaw.org/isotopic-abundances.htm'
    
    
    def __init__(self):
        
        self.setup()
    
    
    def reload(self):
        
        modname = self.__class__.__module__
        mod = __import__(modname, fromlist=[modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)
    
    
    def setup(self):
        """
        Populates the mass database.
        """
        
        self.load_mass_monoiso()
        self.load_freq_iso()
        self.setup_mass_first_iso()
        # self.get_weights() # this does not work at the moment
        self.setup_isotopes()
    
    
    @staticmethod
    def load_masses(url):
        """
        Downloads an HTML table from CIAAW webpage
        and extracts the atomic mass or weight information.

        Parameters
        ----------
        url : str
            URL for the table an HTML table of masses on CIAAW webpage.

        Returns
        -------
        Dict of masses or weights.
        """
        
        c = _curl.Curl(url, silent = False)
        req_masses = c.result
        
        with warnings.catch_warnings():
            # there is a deprecated call in lxml
            warnings.simplefilter('ignore', DeprecationWarning)
            soup_masses = bs4.BeautifulSoup(req_masses, 'lxml')

        masses = {}
        symbol = None
        a = None

        for tr in soup_masses.find_all('tr'):
            
            tr = [td for td in tr.find_all('td')]
            
            if not len(tr):
                
                continue
                
            elif len(tr) == 5:
                
                symbol = tr[1].text.strip()
                masses[symbol] = {}
            
            a = int(_re_nondigit.sub('', tr[-2].text.strip()))
            m = [
                float(_re_nondigit.sub('', i))
                for i in tr[-1].text.split(',')
            ]
            m = sum(m) / len(m)
            masses[symbol][a] = m
        
        masses['proton']   = 1.00727646677
        masses['electron'] = 0.00054857990924
        masses['neutron']  = 1.00866491588
        
        return masses
    

    def load_mass_monoiso(self):
        """
        Obtains monoisotopic masses from CIAAW webpage.
        Stores the result in ``massMonoIso`` module level variable.
        """
        
        self.mass_monoiso = self.load_masses(self.url_masses)
    
    
    def load_freq_iso(self):
        """
        Obtains isotope abundances from CIAAW webpage.
        Stores the result in :py:attr:`.freqIso` attribute of the module.
        """
        
        c = _curl.Curl(self.url_abundances, silent = False)
        req_abundances = c.result.split('\n')
        
        # fixing erroneous HTML from CIAAW:
        for i, l in enumerate(req_abundances[:-1]):
            
            l = l.strip()
            
            if l[-5:] == '</tr>' and req_abundances[i + 1][:3] == '<td':
                
                req_abundances[i + 1] = '<tr>%s' % req_abundances[i + 1]
        
        with warnings.catch_warnings():
            # there is a deprecated call in lxml
            warnings.simplefilter('ignore', DeprecationWarning)
            soup_abundances = bs4.BeautifulSoup(
                '\n'.join(req_abundances),
                'lxml'
            )
        
        freq_iso = {}
        symbol = None
        a = None
        
        for tr in soup_abundances.find_all('tr'):
            
            tr = [td for td in tr.find_all('td')]
            
            if len(tr) == 6:
                
                symbol = tr[1].text.strip()
                freq_iso[symbol] = {}
            
            ai = -3 if len(tr) == 6 else -2
            
            try:
                
                a = int(tr[ai].text.strip())
                p = [
                    float(_re_nondigit.sub('', i))
                    for i in tr[ai + 1].text.split(',')
                ]
                p = sum(p) / len(p)
                freq_iso[symbol][a] = p
                
            except (ValueError, IndexError, KeyError):
                
                continue
        
        self.freq_iso = freq_iso


    def setup_mass_first_iso(self):
        """
        Obtains the masses of the most abundant isotope for each element.
        The result stored in the :py:attr:``mass_first_iso`` attribute.
        """
        
        first_iso = {}
        
        for symbol, isos in iteritems(self.mass_monoiso):
            
            if symbol in self.freq_iso:
                
                try:
                    
                    first_iso[symbol] = (
                        isos[
                            max(
                                iteritems(self.freq_iso[symbol]),
                                key = lambda i: i[1]
                            )[0]
                        ]
                    )
                
                except:
                    
                    continue
        
        self.mass_first_iso = first_iso


    def get_weights(self):
        """
        Obtains atomic weights from CIAAW webpage.
        """
        
        self.weights = self.load_masses(self.url_weights)
    
    
    def setup_isotopes(self):
        """
        Builds a dict of isotopes with their mass and abundance.
        Result stored in :py:attr:``isotopes`` attribute.
        """
        
        isotopes = {}
        
        for elem, isos in iteritems(self.freq_iso):
            
            isotopes[elem] = {}
            
            min_nominal = min(isos.keys(), default = 0)
            
            for nominal, freq in iteritems(isos):
                
                isotopes[elem][nominal - min_nominal] = (
                    self.mass_monoiso[elem][nominal],
                    freq
                )
        
        self.isotopes = isotopes
    
    
    def first_isotope_mass(self, elem):
        """
        Returns exact mass of the highest abundant isotope of an element.

        Parameters
        ----------
        elem : str
            Element.

        Returns
        -------
        Monoisotopic mass of the most abundant isotope of the element.
        """
        
        return (
            self.mass_first_iso[elem]
                if elem in self.mass_first_iso else
            None
        )


    def get_mass(self, elem):
        """
        Returns exact mass of the highest abundant isotope of an element.

        Parameters
        ----------
        elem : str
            Element.

        Returns
        -------
        Monoisotopic mass of the most abundant isotope of the element.
        """
        
        return self.first_isotope_mass(elem)


    def isotope_mass(self, elem, iso):
        """
        Parameters
        ----------
        elem : str
            Element.
        iso : int
            Nominal mass of the isotope.

        Returns
        -------
        Mass of the isotope, None if element or isotope not in database.
        """
        
        return (
            self.mass_monoiso[elem][iso]
                if (
                    elem in self.mass_monoiso and
                    iso in self.mass_monoiso[elem]
                ) else
            None
        )


    def get_weight(self, elem):
        """
        Parameters
        ----------
        elem : str
            Element.

        Returns
        -------
        Atomic weight as float, None if element not in database.
        """
        
        return self.weights[elem] if elem in self.weights else None


def init_db():
    
    globals()['db'] = MassDatabase()


# databases set up at module loading
init_db()


class MassBase(object):
    """
    Represents a mass as a floating point number optionally with a chemical
    formula.
    """
    
    def __init__(
            self,
            formula_mass = None,
            charge = 0,
            isotope = 0,
            **kwargs
        ):
        """
        This class is very similar to ``Formula``. Actually it can be
        initialized either with providing a formula or a mass or
        even element counts as keyword arguments.
        The key difference compared to ``Formula`` is that it behaves
        as a ``float`` i.e. indeed represents a molecular mass, while
        ``Formula`` behaves as a chemical formula i.e. representing
        the counts of elements. If you add two `MassBase` instances
        (or a float) you get a ``float`` while if you add two
        ``Formula`` instances (or a string) you get a new ``Formula``.
        Finally ``Mass`` is able to provide both behaviours but
        adding two ``Mass`` instances will result a new ``Mass``.
        
        Parameters
        ----------
        formula_mass : str,float,NoneType
            Either a string expressing a chemical formula (e.g. H2O) or
            a molecular mass (e.g. 237.1567) or `None` if you provide the
            formula as keyword arguments.
        **kwargs :
            Elements & counts, e.g. ``c = 6, h = 12, o = 6``.
        
        References
        ----------
        Atomic and isotopic mass data obtained from the Commission of
        Isotopic Abundances and Atomic Weights (CIAAW), http://www.ciaaw.org/
        """
        
        self.exmass = db.mass_first_iso
        self.charge = charge
        self.isotope = isotope
        
        if formula_mass is None:
            
            self.formula_from_dict(kwargs)
            
        elif hasattr(formula_mass, 'lower'):
            
            self.formula = formula_mass
            
        elif isinstance(formula_mass, MassBase):
            
            if hasattr(formula_mass, 'mass'):
                
                self.mass = formula_mass.mass
                
            if hasattr(formula_mass, 'formula'):
                
                self.formula = formula_mass.formula
            
            self.mass_calculated = formula_mass.mass_calculated
            
        else:
            
            self.formula = None
        
        if type(formula_mass) is float:
            self.mass = formula_mass
        
        self.calc_mass()
    
    
    def __neg__(self):
        
        return -1 * self.mass
    
    
    def __add__(self, other):
        
        return float(other) + self.mass
    
    
    def __radd__(self, other):
        
        return self.__add__(other)
    
    
    def __iadd__(self, other):
        
        self.mass += float(other)
    
    
    def __sub__(self, other):
        
        return self.mass - float(other)
    
    
    def __rsub__(self, other):
        
        return float(other) - self.mass
    
    
    def __isub__(self, other):
        
        self.mass += float(other)
    
    
    def __truediv__(self, other):
        
        return self.mass / float(other)
    
    
    def __rtruediv__(self, other):
        
        return float(other) / self.mass
    
    
    def __itruediv__(self, other):
        
        self.mass /= float(other)
    
    
    def __mul__(self, other):
        
        return self.mass * float(other)
    
    
    def __rmul__(self, other):
        
        return self.__mul__(other)
    
    
    def __imul__(self, other):
        
        self.mass *= float(other)
    
    
    def __float__(self):
        
        return self.mass
    
    
    def __eq__(self, other):
        
        return abs(self.mass - float(other)) <= 0.01
    
    
    def calc_mass(self):
        """
        Calculates the mass from the formula.
        """
        
        if self.has_formula():
            
            if self.formula == '':
                
                self.mass = 0.0
                self.mass_calculated = True
                
            else:
                
                atoms = (
                    _re_form.findall(self.formula)
                    if not hasattr(self, '_atoms')
                    else iteritems(self._atoms)
                )
                m = 0.0
                for element, count in atoms:
                    count = int(count or '1')
                    m += self.exmass[element] * count
                
                if self.isotope:
                    
                    oms_formula = oms.EmpiricalFormula(self.formula)
                    iso_pattern_gen = oms.CoarseIsotopePatternGenerator(
                        self.isotope + 1
                    )
                    isotopes = oms_formula.getIsotopeDistribution(
                        iso_pattern_gen
                    )
                    this_isotope = list(isotopes.getContainer())[-1].getMZ()
                    m = this_isotope
                    
                m -= self.charge * electron
                self.mass = m
                
                self.mass_calculated = self.has_mass()
            
        else:
            
            self.mass_calculated = False
    
    
    def has_mass(self):
        
        return self.mass > 0.0 or (self.formula == '' and self.mass == 0.0)
    
    
    def has_formula(self):
        
        return self.formula is not None
    
    
    def formula_from_dict(self, atoms):
        """
        Creates a formula from a dict of atom counts.
        
        Parameters
        ----------
        atoms : dict
            Dict of atoms i.e. elements as keys and counts as values.
        """
        
        self.formula = ''.join('%s%u'%(elem.capitalize(), num) \
            for elem, num in iteritems(atoms))
    
    def update_atoms(self):
        """
        Sets the atoms dict by processing the ``formula`` attribute.
        """
        
        if self.formula:
            
            self._atoms = formula_to_atoms(self.formula)
    
    @property
    def atoms(self):
        
        self.update_atoms()
        
        return self._atoms if hasattr(self, '_atoms') else {}
    
    def reload(self):
        
        modname = self.__class__.__module__
        mod = __import__(modname, fromlist=[modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)


class ElementComposition(object):
    """
    Represents the element composition of a certain compound or a group of
    compounds. Can be used to calculate isotopic masses.
    """
    
    def __init__(**kwargs):
        """
        Parameters
        ----------
        **kwargs :
            Elements and their relative proportions.
            E.g. ``C = 1, H = 2.05`` or ``c = 10, h = 22, o = 1``.
        """
        
        self.composition = {}
        
        total = sum(v for v in kwargs.values())
        
        for elem, proportion in iteritems(kwargs):
            
            self.composition[elem.capitalize()] = proportion / total
    
    
    def __getattr__(self, elem):
        
        elem = elem.capitalize()
        
        if elem in self.composition:
            
            return composition[elem]
    
    
    def __getitem__(self, elem):
        
        return self.__getattr__(elem)
    
    
    def __iter__(self):
        
        return iteritems(self.composition)


parts = {
    'water': 'H2O',
    'twowater': 'H4O2',
    'carboxyl': 'COOH',
    'aldehyde': 'CHO',
}


for name, form in parts.items():
    
    setattr(sys.modules[__name__], name, MassBase(form))


def calculate(formula):
    """
    Evaluates a string as formula.
    
    Parameters
    ----------
    formula : str
        Expression as a string e.g. ``HCH3CHOHCOOH - water + electron``.
    
    Returns
    -------
    Mass as float.
    """
    
    result = None
    op = '__add__'
    
    for step in replmi.split(formula):
        
        if step == '-':
            
            op = '__sub__'
            continue
        
        if step == '+':
            
            op = '__add__'
            continue
        
        step = step.strip()
        
        if refloat.match(step):
            step = float(step)
        
        if (
            step in globals() and
            isinstance(globals()[step], (float, int, MassBase))
        ):
            step = globals()[step]
        
        if op is not None:
            
            result = getattr(MassBase(result), op)(MassBase(step))
        
        op = None
    
    return result


#: Synonym of :py:func:`calculate`.
expr = calculate


weight_builtin = {
    "proton": 1.00727646677,
    "electron": 0.00054857990924,
    "neutron": 1.00866491588,
    "H":  1.007825,
    "He": 4.002602,
    "Li": 6.941,
    "Be": 9.012182,
    "B":  10.811,
    "C":  12.0107,
    "N":  14.00674,
    "O":  15.9994,
    "F":  18.9984032,
    "Ne": 20.1797,
    "Na": 22.989768,
    "Mg": 24.3050,
    "Al": 26.981539,
    "Si": 28.0855,
    "P":  30.973762,
    "S":  32.066,
    "Cl": 35.4527,
    "Ar": 39.948,
    "K":  39.0983,
    "Ca": 40.078,
    "Sc": 44.955910,
    "Ti": 47.88,
    "V":  50.9415,
    "Cr": 51.9961,
    "Mn": 54.93805,
    "Fe": 55.847,
    "Co": 58.93320,
    "Ni": 58.6934,
    "Cu": 63.546,
    "Zn": 65.39,
    "Ga": 69.723,
    "Ge": 72.61,
    "As": 74.92159,
    "Se": 78.96,
    "Br": 79.904,
    "Kr": 83.80,
    "Rb": 85.4678,
    "Sr": 87.62,
    "Y":  88.90585,
    "Zr": 91.224,
    "Nb": 92.90638,
    "Mo": 95.94,
    "Tc": 98.0,
    "Ru": 101.07,
    "Rh": 102.90550,
    "Pd": 106.42,
    "Ag": 107.8682,
    "Cd": 112.411,
    "In": 114.82,
    "Sn": 118.710,
    "Sb": 121.757,
    "Te": 127.60,
    "I":  126.90447,
    "Xe": 131.29,
    "Cs": 132.90543,
    "Ba": 137.327,
    "La": 138.9055,
    "Ce": 140.115,
    "Pr": 140.90765,
    "Nd": 144.24,
    "Pm": 145.0,
    "Sm": 150.36,
    "Eu": 151.965,
    "Gd": 157.25,
    "Tb": 158.92534,
    "Dy": 162.50,
    "Ho": 164.93032,
    "Er": 167.26,
    "Tm": 168.93421,
    "Yb": 173.04,
    "Lu": 174.967,
    "Hf": 178.49,
    "Ta": 180.9479,
    "W":  183.85,
    "Re": 186.207,
    "Os": 190.2,
    "Ir": 192.22,
    "Pt": 195.08,
    "Au": 196.96654,
    "Hg": 200.59,
    "Tl": 204.3833,
    "Pb": 207.2,
    "Bi": 208.98037,
    "Po": 209,
    "At": 210,
    "Rn": 222,
    "Fr": 223,
    "Ra": 226.0254,
    "Ac": 227,
    "Th": 232.0381,
    "Pa": 213.0359,
    "U":  238.0289,
    "Np": 237.0482,
    "Pu": 244,
    "Am": 243,
    "Cm": 247,
    "Bk": 247,
    "Cf": 251,
    "Es": 252,
    "Fm": 257,
    "Md": 258,
    "No": 259,
    "Lr": 260,
    "Rf": 261,
    "Db": 262,
    "Sg": 263,
    "Bh": 262,
    "Hs": 265,
    "Mt": 266,
}

isotopes_builtin = {
    "H2": 2.01410178,
    "H3": 3.0160492,
    "C13": 13.003355,
    "N15": 15.000109,
    "O17": 16.999132,
    "O18": 17.999160,
    "S33": 32.971458,
    "S34": 33.967867,
    "S35": 35.967081
}

iso_freq_builtin = {
    "H2": 0.000115,
    "H3": 0.0,
    "C13": 0.0107,
    "N15": 0.0068,
    "O17": 0.00038,
    "O18": 0.00205,
    "S33": 0.0076,
    "S34": 0.0429,
    "S35": 0.0002
}
