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

"""
Classes representing fragment ions in MS2 spectra.
"""

from __future__ import print_function
from future.utils import iteritems
from past.builtins import xrange, range, reduce

import re
import sys
import imp

import lipyd.mass as mass
import lipyd.metabolite as metabolite


fattyfragments = set([])


class AdductCalculator(object):
    
    def __init__(self):
        self.reform = re.compile(r'([A-Za-z][a-z]*)([0-9]*)')
        self.init_counts()
    
    def init_counts(self):
        self.counts = self.counts if hasattr(self, 'counts') else {}
    
    def formula2counts(self, formula):
        for elem, num in self.reform.findall(formula):
            yield elem, int(num or '1')
    
    def remove(self, formula):
        for elem, num in self.formula2counts(formula):
            self.add_atoms(elem, -num)
    
    def add(self, formula):
        for elem, num in self.formula2counts(formula):
            self.add_atoms(elem, num)
    
    def add_atoms(self, elem, num):
        self.counts[elem] = num if elem not in self.counts \
            else self.counts[elem] + num


class FattyFragment(metabolite.AbstractSubstituent):
    
    def __init__(
            self,
            head = '',
            minus = '',
            charge = 0,
            name = 'UnknownFattyFragment%s',
            attrs = None,
            mode = None,
            **kwargs
        ):
        
        self.mode = mode
        
        def getname(parent):
            
            return '[%s]%s' % (
                parent.name % ('(C%u:%u)' % (parent.c, parent.u)),
                'NL' if parent.name[:2] == 'NL' else
                '+'  if parent.charge == 1 else
                '-'  if parent.charge == -1 else
                ''
            )
        
        cminus = dict(
            (elem, -cnt) for elem, cnt in
            iteritems(metabolite.formula.formula2atoms(minus))
        )
        
        cminus['H'] = cminus['H'] - 2 if 'H' in cminus else -2
        
        attrs = attrs or {}
        attrs['typ'] = (
            attrs['typ']
            if 'typ' in attrs else
            lambda parent: (
                parent.name % ''
                if '%s' in parent.name else
                parent.name
            )
        )
        
        metabolite.AbstractSubstituent.__init__(
            self,
            cores = [head],
            counts = cminus,
            names = name,
            charges = charge,
            getname = getname,
            # we keep this 0 to avoid complicating further
            valence = 0,
            attrs = attrs,
            **kwargs
        )
    
    def iterfraglines(self):
        
        for fr in self:
            
            yield [
                fr.mass, fr.name, fr.attrs.typ,
                self.c, self.u, self.charge
            ]

class FattyFragmentOld(mass.MassBase, AdductCalculator):
    
    def __init__(self, charge, c = 3, unsat = 0,
        minus = None, plus = None, isotope = 0, name = None, hg = None):
        self.c = c
        self.unsat = unsat
        self.minus = minus or []
        self.plus  = plus or []
        self.hg = hg or []
        self.init_counts()
        self.add_atoms('C', self.c)
        self.add_atoms('H', self.c * 2)
        self.add_atoms('O', 2)
        self.add_atoms('H', self.unsat * - 2 )
        AdductCalculator.__init__(self)
        for formula in self.minus:
            self.remove(formula)
        for formula in self.plus:
            self.add(formula)
        mass.MassBase.__init__(self, charge = charge,
            isotope = isotope, **self.counts)
        self.set_name(name)
    
    def set_name(self, name):
        if name is None: name = ''
        self.name = '[%s(C%u:%u)%s%s]%s' % (
            name,
            self.c,
            self.unsat,
            self.adduct_str(),
            '' if self.isotope == 0 else 'i%u' % self.isotope,
            self.charge_str()
        )
    
    def charge_str(self):
        return '%s%s' % (
                ('%u' % abs(self.charge)) if abs(self.charge) > 1 else '',
                '-' if self.charge < 0 else '+'
            ) if self.charge != 0 else 'NL'
    
    def adduct_str(self):
        return '%s%s' % (
                ('-' if self.charge < 0 else '') \
                    if not self.minus  else \
                        '-%s' % ('-'.join(self.minus)),
                ('+' if self.charge > 0 else '') \
                    if not self.plus  else \
                        '+%s' % ('+'.join(self.plus))
            )
    
    def get_fragline(self):
        return [self.mass, self.name,
            '[M%s]%s' % (self.adduct_str(), self.charge_str()),
            ';'.join(self.hg)]


class FattyFragmentFactory(object):
    
    docs = {
        'LysoPEAlkyl':
            """
            from massbank.jp:
            [lyso PE(alkenyl-18:0,-)]- 464.3140997565 -417 C23H47NO6P-
            """,
        'LysoPE':
            """
            from massbank.jp:
            [lyso PE(18:0,-)]- 480.3090143786 -476 C23H47NO7P-
            """,
        'LysoPCAlkyl':
            """
            from massbank.jp:
            [lyso PC(alkenyl-18:0,-)]- 492.3453998849 -436 C25H51NO6P-
            
            from massbank.jp:
            [lyso PC(alkyl-18:0,-)]- 494.3610499491 -143 C25H53NO6P-
            """,
        'LysoPC':
            """
            from massbank.jp:
            [lyso PC(18:0,-)]- 508.340314507 -373 C25H51NO7P-
            """,
        'LysoPA_mH2O':
            """
            from massbank.eu:
            https://massbank.eu/MassBank/jsp/RecordDisplay.jsp?id=UT002963&dsn=CHUBU
            [lyso PS(18:0,-)-H2O]- 419.2562505266 -610 C21H40O6P-
            
            Frega 2007 Fig 15a
            """,
        'LysoPI':
            """
            from massbank:
            [lyso PI(-,18:0)]- 599.3196386444 -432 C27H52O12P-
            """,
        'LysoPIAlkyl':
            """
            from massbank.jp:
            [lyso PI(alkyl-16:1,-)-H2O]- 537.2828592076 -228 C25H46O10P-
            
            from this, derived 18:0-:
            [lyso PI(alkyl-18:0,-)]- 585.3403740858 -228 C27H54O11P-
            """,
        'LysoPG':
            """
            from massbank:
            [lyso PG(18:0,-)]- 511.3035946497 -495 C24H48O9P-
            """,
        'LysoPGAlkyl':
            """
            from massbank:
            [lyso PG(18:0,-)]- 511.3035946497 -495 C24H48O9P-
            """,
        'LysoPAAlkyl':
            """
            from Characterization of Phospholipid
            Molecular Species by Means of HPLC-Tandem Mass Spectrometry:
            [lyso PA(18:1,-)]- 435.2 C21H40O7P-
            """,
        'LysoPA':
            """
            from Characterization of Phospholipid
            Molecular Species by Means of HPLC-Tandem Mass Spectrometry:
            [lyso PA(18:1,-)]- 435.2 C21H40O7P-
            
            from massbank.jp:
            https://massbank.eu/MassBank/jsp/RecordDisplay.jsp?id=UT002963&dsn=CHUBU
            [lyso PS(18:0,-)]- 437.2668152129 -358 C21H42O7P-
            
            Frega 2007 Fig 15a
            """,
        'FA_mO_pC2H2NH2':
            """
            https://metlin.scripps.edu/metabo_info.php?molid=6214
            [Cer-FA(C8:0)]- 168.1382904 C8H14O1N1-
            """,
        'FA_mH2O_mH':
            """
            209 at C14:0 FA, 263 at C18:1
            """,
        'FA_mO_pNH2':
            """
            226 at C14:0 FA, 280 at C18:1
            """,
        'FA_mH':
            """
            227 at C14:0 FA, 281 at C18:1
            """,
        'FAL_mH':
            """
            18:1 = 267.2693393
            """
    }
    
    # class name: (head formula, minus, charge, name, headgroups, type)
    # TODO make this an input file instead of dict
    # in comments masses at 18:0
    param_neg = {
        'LysoPE':           ('C5H9O7NH2P', '', -1, 'LysoPE%s', ['PE'], 'fa'),
        'LysoPEAlkyl':      ('C5H11O6NH2P', '', -1, 'LysoPEAlkyl%s', ['PE'],
                             'fal'),
        'LysoPCAlkyl':      ('C7H17O6NP', '', -1, 'LysoPCAlkyl%s', ['PC'],
                             'fal'),
        'LysoPC':           ('C7H15O7NP', '', -1, 'LysoPC%s', ['PC'], 'fa'),
        'LysoPI':           ('C9H16O12P', '', -1, 'LysoPI%s', ['PI'], 'fa'),
        'LysoPIAlkyl':      ('C9H18O11P', '', -1, 'LysoPIAlkyl%s', ['PI'],
                             'fal'),
        'LysoPG':           ('C6H12O9P', '', -1, 'LysoPG%s', ['PG'], 'fa'),
        'LysoPGAlkyl':      ('C6H14O8P', '', -1, 'LysoPGAlkyl%s', ['PG'],
                             'fal'),
        'LysoPA':           ('C3H6O7P',  '', -1, 'LysoPA%s', ['PA', 'PS'],
                             'fa'),
        'LysoPAAlkyl':      ('C3H8O6P',  '', -1, 'LysoPAAlkyl%s',
                             ['PA', 'PS'], 'fal'),
        'LysoPA_mH2O':      ('C3H4O6P',  '', -1, 'LysoPA%s-H2O',
                             ['PA', 'PS'], 'fa'),
        # former CerFA
        'FA_mO_pC2H2NH2':   ('C2H2ON', '',   -1, 'FA%s-O+C2H2NH2', ['Cer'],
                             'fa'), # 308.2959
        # former CerFAminusC2H5N
        'FA_mH2O_mH':       ('O', 'H3', -1, 'FA%s-H2O-H', ['Cer'],
                             'fa'), # 265.2537
        # former CerFAminusC
        'FA_mO_pNH2':       ('NO', '', -1, 'FA%s-O+NH2', ['Cer'], 'fa'),
        # former CerSphiMinusN
        'Sph_mC2H4_mNH2_mH2O':
                            ('O', 'C2H5', -1, 'Sph%s-C2H4-NH2-H2O', ['Cer'],
                             'sphingo'), # 239.2380
        # former CerSphiMinusNO
        'Sph_mH2O_mNH2_m2H':
                            ('O', 'H3', -1, 'Sph%s-H2O-NH2-2H', ['Cer'],
                             'sphingo'), # 265.2537
        # former CerSphi
        'Sph_mC2H4_m3H':    ('NO2', 'C2H4', -1, 'Sph%s-C2H4-3H', ['Cer'],
                             'sphingo'), # 270.2436
        # former CerFAminusN, FAminusH
        'FA_mH':            ('O2', 'H', -1, 'FA%s-H', [], 'fa'), # 283.2643
        # former FAAlkylminusH
        'FAL_mH':           ('OH', '', -1, 'FAL%s-H', [], 'fal'),
        # fatty acyl fragments for hydroxyacyl ceramides:
        'FA_pC2H3_pNH2_pH2O': ('C2NH4O3', '', -1, 'FA%s+C2H3+NH2+H2O', [],
                                 'fa'), # 342.3014
        'FA_pC2H3_pNH2_pO': ('C2NH2O3', '', -1, 'FA%s+C2H3+NH2+O', [],
                              'fa'), # 340.2857
        'FA_pC2H3_pNH2':    ('C2NH2O2', '', -1, 'FA%s+C2H3+NH2', [],
                             'fa'), # 324.2908
        'FA_pO':            ('O3', 'H', -1, 'FA%s+O', [], 'fa'), # 299.2592
        'FA_pC2H3_pNH2_m2H':
                            ('C2NO2', '', -1, 'FA%s+C2H3+NH2-2H', [],
                             'fa'), # 322.2752
        'FA_pC2H3_pNH2_mH2O': ('C2NO', '', -1, 'FA%s+C2H3+NH2-H2O', [],
                               'fa'), # 306.2802
        'FA_pNH2_m2H':      ('NO2', '', -1, 'FA%s+NH2-2H', [],
                             'fa'), # 298.2752
        'FA_m3H':           ('O2', 'H3', -1, 'FA%s-3H', [], 'fa'), # 281.2486
        'FA_mCO_m3H':       ('O', 'CH3', -1, 'FA%s-CO-3H', [],
                             'fa'), # 253.2537
        # fatty acyl fragments from ceramides
        'FA_pC2H3_pNH2_p2H': ('C2NH4O2', '', -1, 'FA%s+C2H3+NH2+2H', [],
                              'fa'), # 326.3065
        'FA_pC2H3_pNH2_mO': ('C2NH2O', '', -1, 'FA%s+C2H3+NH2-O', [],
                               'fa'), # 308.2959
        'FA_pNH2_mH2O_mH':  ('ONH', '', -1, 'FA%s+NH2-H2O-H', [],
                             'fa'), # 282.2802
        'FA_pC3H7_mNH2':    ('O2NC3H4', '', -1, 'FA%s+C3H7+NH2', [],
                             'fa'), # 338.3064
        # sphingosine long chain base fragments from t-ceramides
        'Sph_pO':           ('NH2O3', '', -1, 'Sph%s+O', ['Cer'],
                             'sphingo'), # 316.2857
        'Sph_mCH2_m3H':     ('NO2', 'CH2', -1, 'Sph%s-CH2-3H', ['Cer'],
                             'sphingo'), # 284.2595
        'Sph_mNH2_m3H':     ('O2', 'H3', -1, 'Sph%s-NH2-3H', ['Cer'],
                             'sphingo'), # 281.2486
        'Sph_mC2H4_mNH2_m2H':
                            ('O2', 'C2H5', -1, 'Sph%s-C2H4-NH2-2H', ['Cer'],
                             'sphingo'), # 255.2330
        'Sph_mCH2_mNH2_m4H':
                            ('O2', 'CH5', -1, 'Sph%s-CH2-NH2-4H', ['Cer'],
                             'sphingo'), # 267.2330
        # sphingosine long chain base fragments from DH-ceramides
        'Sph_mH_N':         ('NH2O2', '', -1, 'Sph%s-H', ['Cer'],
                             'sphingo'), # 300.2908
        'Sph_mCH2_mH2O_mH': ('NO', 'CH2', -1, 'Sph%s-CH2-H2O-H', ['Cer'],
                             'sphingo'), # 268.2646
        # sphingosine long chain base fragments from d-ceramides
        'Sph_mCH2_mOH':     ('NO', 'C', -1, 'Sph%s-CH2-OH', ['Cer'],
                             'sphingo'), # 270.2802
        # sphingosine long chain base fragments from sphingomyelin
        'Sph_pPCh_mCH3':    ('N2H12O5PC4', '', -1, 'Sph%s+PCh-CH3', ['SM'],
                             'sphingo'), # 449.3150
        # sphingosine fragments from sphingosine
        'Sph_mC2H3_mNH2_mOH':
                            ('O', 'C2H3', -1, 'Sph%s-C2H3-NH2-OH', ['Sph'],
                             'sphingo'), # 241.2537
        # fatty acyl fragments from ceramide-1-phosphate
        'FA_mH2O':          ('O', 'H2', -1, 'FA%s-H2O', ['CerP'],
                             'fa'), # 266.2615
        'FA':               ('O2', '', -1, 'FA%s-', ['CerP'],
                             'fa'), # 284.2721
    }
    
    param_pos = {
        ### neutral losses
        'NLFA':             ('O2', '', 0, 'NL FA%s', [], 'fa'),
        # NLFAminusH2O
        'NLFA_mH2O':        ('O', 'H2', 0, 'NL FA%s-H20', [], 'fa'),
        # NLFAplusOH
        'NLFA_pOH':         ('O3H', '', 0, 'NL FA%s+OH', [], 'fa'),
        # NLFAplusNH3
        'NLFA_pNH3':        ('O2NH3', '', 0, 'NL FA%s+NH3', [], 'fa'),
        ### positive
        # FAminusO
        'FA_mOH':           ('O', 'H', 1, 'FA%s-OH', [], 'fa'),
        # FAplusGlycerol
        'FA_pGlycerol_mOH': ('C3H5O3', '', 1, 'FA%s+Glycerol-OH',
                             ['PG', 'BMP', 'DAG', 'LysoPE', 'LysoPC'], 'fa'),
        # SphingosineBase
        'Sph_pH':           ('NH4O2', '', 1, 'Sph%s+H',
                             ['Sph', 'SM', 'Cer', 'HexCer'],
                             'sphingo'), # 302.3053
        # sphingosine long chain base fragments from ceramide-1-phosphate
        'Sph_m2xH2O':       ('N', '', 1, 'Sph%s-2xH2O',
                             ['CerP', 'SphP', 'SM', 'MSph', 'M2Sph', 'dCer',
                              'DHCer'],
                             'sphingo'), # 266.2842
        # fatty acyl fragments from ceramide-1-phosphate
        'FA_pNH_pC2H2_mOH': ('ONC2H2', '', 1, 'FA%s+NH+C2H2-OH', [],
                             'fa'), # 308.2948
        'FA_pNH2_mO':       ('ONH2', '', 1, 'FA%s+NH2-O', ['dCer', 'DHCer'],
                             'fa'), # 284.2948
        # sphingosine long chain base fragments from hexosyl-t-ceramide
        'Sph_mH_P':         ('NH2O2', '', 1, 'Sph%s-H', ['HexCer', 'tCer'],
                             'sphingo'), # 300.2897
        'Sph_pH2O_mH':      ('NH4O3', '', 1, 'Sph%s+H2O-H', ['HexCer'],
                             'sphingo'), # 318.3003
        'Sph_mH2O_mH':      ('NO', '', 1, 'Sph%s-H2O-H', ['HexCer', 'tCer'],
                             'sphingo'), # 282.2791
        'Sph_mH2O_pH':      ('NOH2', '', 1, 'Sph%s-H2O+H',
                             ['SphP', 'MSph', 'M2Sph', 'DHCer'],
                             'sphingo'), # 284.2948
        # sphingosine fragments from N-(n)methyl-safingols
        'Sph_mH2O_pCH3':    ('NOCH4', '', 1, 'Sph%s-H2O+CH3',
                             ['MSph', 'M2Sph'], 'sphingo'), # 298.3104
        'Sph_mH2O_p2xCH3_mH':
                            ('NOC2H6', '', 1, 'Sph%s-H2O+2xCH3-H',
                             ['M2Sph'], 'sphingo'), # 312.3261
        # sphinosine fragments from keto-sphingosine
        'Sph_mC_mH2O_mH':   ('NO', 'C', 1, 'Sph%s-C-H2O-H', ['KSph'],
                             'sphingo'), # 270.2791
        'Sph_mC_mO_mH2O_mH':
                            ('N', 'C', 1, 'Sph%s-C-O-H2O-H', ['KSph', 'dCer'],
                             'sphingo'), # 254.2842
        'Sph_mNH2_mH2O_m2H':
                            ('O', 'H3', 1, 'Sph%s-C-H2O-2H', ['KSph'],
                             'sphingo'), # 265.2526
        'Sph_mC_m2xH2O':    ('N', 'CH2', 1, 'Sph%s-C-2xH2O', ['KSph'],
                             'sphingo'), # 252.2686
        # sphingosine fragments from d-ceramide
        'Sph_mO_pH':        ('NOH4', '', 1, 'Sph%s-O+H',
                             ['dCer'], 'sphingo'), # 286.3104
    }
    
    def __init__(self):
        
        mod = sys.modules[__name__]
        
        for mode in ('neg', 'pos'):
            
            param = getattr(self, 'param_%s' % mode)
            
            for name, par in iteritems(param):
                
                exec(
                    (
                        'def __init__(self, **kwargs):\n'
                        '    FattyFragment.__init__(\n'
                        '        self,\n'
                        '        head = \'%s\',\n'
                        '        minus = \'%s\',\n'
                        '        name = \'%s\',\n'
                        '        charge = %u,\n'
                        '        mode = \'%s\',\n'
                        '        attrs = {\n'
                        '            \'headgroups\': %s,\n'
                        '            \'type\': \'%s\'\n'
                        '        },\n'
                        '        **kwargs\n'
                        '    )\n'
                    ) % (
                        par[0],
                        par[1],
                        par[3],
                        par[2],
                        mode,
                        str(par[4]),
                        par[5]
                    ),
                    mod.__dict__,
                    mod.__dict__
                )
                
                if name in self.docs:
                    
                    mod.__dict__['__init__'].__doc__ = self.docs[name]
                
                cls = type(
                    name,
                    (FattyFragment, ),
                    {'__init__': mod.__dict__['__init__']}
                )
                
                # these we set as attributes of the class
                # in order to be able later to inspect them
                # and get information about the class
                for attr, val in (
                    ('typ', par[5]),
                    ('headgroups', par[4]),
                    ('mode', mode)
                ):
                    
                    setattr(cls, attr, val)
                
                setattr(mod, name, cls)
                
                fattyfragments.add(name)
        
        delattr(mod, '__init__')


_factory = FattyFragmentFactory()
del _factory


class FAFragSeries(object):
    
    def __init__(self, typ, charge, cmin = 2, unsatmin = 0,
        cmax = 36, unsatmax = 6,
        minus = [], plus = [], **kwargs):
        
        for attr, val in iteritems(locals()):
            setattr(self, attr, val)
        
        self.fragments = []
        if self.unsatmax is None: self.unsatmax = self.unsatmin
        for unsat in xrange(self.unsatmin, self.unsatmax + 1):
            
            this_cmin = max(self.cmin, unsat * 2 + 1)
            if this_cmin <= cmax:
                
                for cnum in xrange(this_cmin, cmax + 1):
                    
                    self.fragments.append(
                        self.typ(charge = charge, c = cnum, unsat = unsat,
                            minus = self.minus, plus = self.plus, **kwargs)
                    )
    
    def __iter__(self):
        for fr in self.fragments:
            yield fr
    
    def itermass(self):
        for fr in self.fragments:
            yield fr.mass
    
    def iterfraglines(self):
        for fr in self.fragments:
            yield fr.get_fragline()
