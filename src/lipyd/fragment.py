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
#

from __future__ import print_function
from future.utils import iteritems
from past.builtins import xrange, range, reduce

import re
import sys
import imp
import collections

import lipyd.mass as mass
import lipyd.metabolite as metabolite


fattyfragments = set([])


ChainFragParam = collections.namedtuple(
    'ChainFragParam',
    ['plus', 'minus', 'charge', 'name', 'chaintype', 'constraints'] )
ChainFragParam.__new__.__defaults__ = ((),)


FragConstraint = collections.namedtuple(
    'FragConstraint',
    ['hg', 'family', 'sub', 'sph', 'oh', 'chaintype']
)
FragConstraint.__new__.__defaults__ = (None, None, None, None, 0, None)


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
            ionmode = None,
            **kwargs
        ):
        
        self.ionmode = ionmode
        
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
        attrs['fragtype'] = (
            attrs['fragtype']
            if 'fragtype' in attrs else
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
                fr.mass, fr.name, fr.attrs.fragtype, fr.attrs.chaintype,
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
            from massbank.jp:3
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
        'LysoPE':
            ChainFragParam(
                plus = 'C5H9O7NH2P',
                minus = '',
                charge = -1,
                name = 'LysoPE%s',
                constraints = (
                    FragConstraint(
                        hg = 'PE',
                        chaintype = 'FA',
                    ),
                ),
                chaintype = 'FA'
            ),
        'LysoPEAlkyl':
            ChainFragParam(
                plus = 'C5H11O6NH2P',
                minus = '',
                charge = -1,
                name = 'LysoPEAlkyl%s',
                constraints = (
                    FragConstraint(
                        hg = 'PE',
                        chaintype = 'FAL',
                    ),
                ),
                chaintype = 'FAL'
            ),
        'LysoPEAlkyl_mH2O':
            ChainFragParam(
                plus = 'C5H9O5NH2P',
                minus = '',
                charge = -1,
                name = 'LysoPEAlkyl%s-H2O',
                constraints = (
                    FragConstraint(
                        hg = 'PE',
                        chaintype = 'FAL',
                    ),
                ),
                chaintype = 'FAL'
            ),
        'LysoPCAlkyl':
            ChainFragParam(
                plus = 'C7H17O6NP',
                minus = '',
                charge = -1,
                name = 'LysoPCAlkyl%s',
                constraints = (
                    FragConstraint(
                        hg = 'PC',
                        chaintype = 'FAL',
                    ),
                ),
                chaintype = 'FAL'
            ),
        'LysoPC':
            ChainFragParam(
                plus = 'C7H15O7NP',
                minus = '',
                charge = -1,
                name = 'LysoPC%s',
                constraints = (
                    FragConstraint(
                        hg = 'PC',
                        chaintype = 'FA',
                    ),
                ),
                chaintype = 'FA'
            ),
        'LysoPI':
            ChainFragParam(
                plus = 'C9H16O12P',
                minus = '',
                charge = -1,
                name = 'LysoPI%s',
                constraints = (
                    FragConstraint(hg = 'PI', chaintype = 'FA'),
                ),
                chaintype = 'FA'
            ),
        'LysoPI_mH2O':
            ChainFragParam(
                plus = 'C9H14O11P',
                minus = '',
                charge = -1,
                name = 'LysoPI%s-H2O',
                constraints = (
                    FragConstraint(hg = 'PI', chaintype = 'FA'),
                ),
                chaintype = 'FA'
            ),
        'LysoPIAlkyl':
            ChainFragParam(
                plus = 'C9H18O11P',
                minus = '',
                charge = -1,
                name = 'LysoPIAlkyl%s',
                constraints = (
                    FragConstraint(hg = 'PI', chaintype = 'FAL'),
                ),
                chaintype = 'FAL'
            ),
        'LysoPG':
            ChainFragParam(
                plus = 'C6H12O9P',
                minus = '',
                charge = -1,
                name = 'LysoPG%s',
                constraints = (
                    FragConstraint(hg = 'PG', chaintype = 'FA'),
                    FragConstraint(hg = 'BMP', chaintype = 'FAL'),
                ),
                chaintype = 'FA'
            ),
        'LysoPG_mH2O':
            ChainFragParam(
                plus = 'C6H10O8P',
                minus = '',
                charge = -1,
                name = 'LysoPG%s-H2O',
                constraints = (
                    FragConstraint(hg = 'PG', chaintype = 'FA'),
                    FragConstraint(hg = 'BMP', chaintype = 'FAL'),
                ),
                chaintype = 'FA'
            ),
        'LysoPGAlkyl':
            ChainFragParam(
                plus = 'C6H14O8P',
                minus = '',
                charge = -1,
                name = 'LysoPGAlkyl%s',
                constraints = (
                    FragConstraint(hg = 'PG', chaintype = 'FAL'),
                ),
                chaintype = 'FAL'
            ),
        'LysoPA':
            ChainFragParam(
                plus = 'C3H6O7P',
                minus = '',
                charge = -1,
                name = 'LysoPA%s',
                constraints = (
                    FragConstraint(hg = 'PA', chaintype = 'FA'),
                    FragConstraint(hg = 'PS', chaintype = 'FA'),
                ),
                chaintype = 'FA'
            ),
        'LysoPAAlkyl':
            ChainFragParam(
                plus = 'C3H8O6P',
                minus = '',
                charge = -1,
                name = 'LysoPAAlkyl%s',
                constraints = (
                    FragConstraint(hg = 'PA', chaintype = 'FAL'),
                    FragConstraint(hg = 'PS', chaintype = 'FAL'),
                ),
                chaintype = 'FAL'
            ),
        'LysoPA_mH2O':
            ChainFragParam(
                plus = 'C3H4O6P',
                minus = '',
                charge = -1,
                name = 'LysoPA%s-H2O',
                constraints = (
                    FragConstraint(hg = 'PA', chaintype = 'FA'),
                    FragConstraint(hg = 'PS', chaintype = 'FA'),
                ),
                chaintype = 'FA'
            ),
        # former CerFA
        'FA_mO_pC2H2NH2':
            ChainFragParam(
                plus = 'C2H2ON',
                minus = '',
                charge = -1,
                name = 'FA%s-O+C2H2NH2',
                constraints = (
                    FragConstraint(hg = 'Cer', chaintype = 'FA'),
                ),
                chaintype = 'FA'
            ), # 308.2959
        # former CerFAminusC2H5N
        'FA_mH2O_mH':
            ChainFragParam(
                plus = 'O',
                minus = 'H3',
                charge = -1,
                name = 'FA%s-H2O-H',
                constraints = (
                    FragConstraint(hg = 'Cer', chaintype = 'FA'),
                ),
                chaintype = 'FA'
            ), # 265.2537
        # former CerFAminusC
        'FA_mO_pNH2':
            ChainFragParam(
                plus = 'NO',
                minus = '',
                charge = -1,
                name = 'FA%s-O+NH2',
                constraints = (
                    FragConstraint(hg = 'Cer', chaintype = 'FA'),
                ),
                chaintype = 'FA'
            ),
        # former CerSphiMinusN
        'Sph_mC2H4_mNH2_mH2O':
            ChainFragParam(
                plus = 'O',
                minus = 'C2H5',
                charge = -1,
                name = 'Sph%s-C2H4-NH2-H2O',
                constraints = (
                    FragConstraint(hg = 'Cer', chaintype = 'Sph'),
                ),
                chaintype = 'Sph'
            ), # 239.2380
        # former CerSphiMinusNO
        'Sph_mH2O_mNH2_m2H':
            ChainFragParam(
                plus = 'O',
                minus = 'H3',
                charge = -1,
                name = 'Sph%s-H2O-NH2-2H',
                constraints = (
                    FragConstraint(hg = 'Cer', chaintype = 'Sph'),
                ),
                chaintype = 'Sph'
            ), # 265.2537
        # former CerSphi
        'Sph_mC2H4_m3H':
            ChainFragParam(
                plus = 'NO2',
                minus = 'C2H4',
                charge = -1,
                name = 'Sph%s-C2H4-3H',
                constraints = (
                    FragConstraint(hg = 'Cer', chaintype = 'Sph'),
                ),
                chaintype = 'Sph'
            ), # 270.2436
        # former CerFAminusN, FAminusH
        'FA_mH':
            ChainFragParam(
                plus = 'O2',
                minus = 'H',
                charge = -1,
                name = 'FA%s-H',
                constraints = (
                    FragConstraint(chaintype = 'FA'),
                ),
                chaintype = 'FA'
            ), # 283.2643
        # former FAAlkylminusH
        'FAL_mH':
            ChainFragParam(
                plus = 'OH',
                minus = '',
                charge = -1,
                name = 'FAL%s-H',
                constraints = (
                    FragConstraint(chaintype = 'FAL'),
                ),
                chaintype = 'FAL'
            ),
        # fatty acyl fragments for hydroxyacyl ceramides:
        'FA_pC2H3_pNH2_pH2O':
            ChainFragParam(
                plus = 'C2NH4O3',
                minus = '',
                charge = -1,
                name = 'FA%s+C2H3+NH2+H2O',
                constraints = (
                    FragConstraint(family = 'SL', chaintype = 'FA'),
                ),
                chaintype = 'FA'
            ), # 342.3014
        'FA_pC2H3_pNH2_pO':
            ChainFragParam(
                plus = 'C2NH2O3',
                minus = '',
                charge = -1,
                name = 'FA%s+C2H3+NH2+O',
                constraints = (
                    FragConstraint(family = 'SL', chaintype = 'FA'),
                ),
                chaintype = 'FA'
            ), # 340.2857
        'FA_pC2H3_pNH2':
            ChainFragParam(
                plus = 'C2NH2O2',
                minus = '',
                charge = -1,
                name = 'FA%s+C2H3+NH2',
                constraints = (
                    FragConstraint(family = 'SL', chaintype = 'FA'),
                ),
                chaintype = 'FA'
            ), # 324.2908
        'FA_pO':
            ChainFragParam(
                plus = 'O3',
                minus = 'H',
                charge = -1,
                name = 'FA%s+O',
                constraints = (
                    FragConstraint(chaintype = 'FA'),
                ),
                chaintype = 'FA'
            ), # 299.2592
        'FA_pC2H3_pNH2_m2H':
            ChainFragParam(
                plus = 'C2NO2',
                minus = '',
                charge = -1,
                name = 'FA%s+C2H3+NH2-2H',
                constraints = (
                    FragConstraint(family = 'SL', chaintype = 'FA'),
                ),
                chaintype = 'FA'
            ), # 322.2752
        'FA_pC2H3_pNH2_mH2O':
            ChainFragParam(
                plus = 'C2NO',
                minus = '',
                charge = -1,
                name = 'FA%s+C2H3+NH2-H2O',
                constraints = (
                    FragConstraint(family = 'SL', chaintype = 'FA'),
                ),
                chaintype = 'FA'
            ), # 306.2802
        'FA_pNH2_m2H':
            ChainFragParam(
                plus = 'NO2',
                minus = '',
                charge = -1,
                name = 'FA%s+NH2-2H',
                constraints = (
                    FragConstraint(family = 'SL', chaintype = 'FA'),
                ),
                chaintype = 'FA'
            ), # 298.2752
        'FA_m3H':
            ChainFragParam(
                plus = 'O2',
                minus = 'H3',
                charge = -1,
                name = 'FA%s-3H',
                constraints = (
                    FragConstraint(chaintype = 'FA'),
                ),
                chaintype = 'FA'
            ), # 281.2486
        'FA_mCO_m3H':
            ChainFragParam(
                plus = 'O',
                minus = 'CH3',
                charge = -1,
                name = 'FA%s-CO-3H',
                constraints = (
                    FragConstraint(chaintype = 'FA'),
                ),
                chaintype = 'FA'
            ), # 253.2537
        # fatty acyl fragments from ceramides
        'FA_pC2H3_pNH2_p2H':
            ChainFragParam(
                plus = 'C2NH4O2',
                minus = '',
                charge = -1,
                name = 'FA%s+C2H3+NH2+2H',
                constraints = (
                    FragConstraint(family = 'SL', chaintype = 'FA'),
                ),
                chaintype = 'FA'
            ), # 326.3065
        'FA_pC2H3_pNH2_mO':
            ChainFragParam(
                plus = 'C2NH2O',
                minus = '',
                charge = -1,
                name = 'FA%s+C2H3+NH2-O',
                constraints = (
                    FragConstraint(family = 'SL', chaintype = 'FA'),
                ),
                chaintype = 'FA'
            ), # 308.2959
        'FA_pNH2_mH2O_mH':
            ChainFragParam(
                plus = 'ONH',
                minus = '',
                charge = -1,
                name = 'FA%s+NH2-H2O-H',
                constraints = (
                    FragConstraint(family = 'SL', chaintype = 'FA'),
                ),
                chaintype = 'FA'
            ), # 282.2802
        'FA_pC3H7_mNH2':
            ChainFragParam(
                plus = 'O2NC3H4',
                minus = '',
                charge = -1,
                name = 'FA%s+C3H7+NH2',
                constraints = (
                    FragConstraint(family = 'SL', chaintype = 'FA'),
                ),
                chaintype = 'FA'
            ), # 338.3064
        # sphingosine long chain base fragments from t-ceramides
        'Sph_pO':
            ChainFragParam(
                plus = 'NH2O3',
                minus = '',
                charge = -1,
                name = 'Sph%s+O',
                constraints = (
                    FragConstraint(hg = 'Cer', chaintype = 'Sph'),
                ),
                chaintype = 'Sph'
            ), # 316.2857
        'Sph_mCH2_m3H':
            ChainFragParam(
                plus = 'NO2',
                minus = 'CH2',
                charge = -1,
                name = 'Sph%s-CH2-3H',
                constraints = (
                    FragConstraint(hg = 'Cer', chaintype = 'Sph'),
                ),
                chaintype = 'Sph'
            ), # 284.2595
        'Sph_mNH2_m3H':
            ChainFragParam(
                plus = 'O2',
                minus = 'H3',
                charge = -1,
                name = 'Sph%s-NH2-3H',
                constraints = (
                    FragConstraint(hg = 'Cer', chaintype = 'Sph'),
                ),
                chaintype = 'Sph'
            ), # 281.2486
        'Sph_mC2H4_mNH2_m2H':
            ChainFragParam(
                plus = 'O2',
                minus = 'C2H5',
                charge = -1,
                name = 'Sph%s-C2H4-NH2-2H',
                constraints = (
                    FragConstraint(hg = 'Cer', chaintype = 'Sph'),
                ),
                chaintype = 'Sph'
            ), # 255.2330
        'Sph_mCH2_mNH2_m4H':
            ChainFragParam(
                plus = 'O2',
                minus = 'CH5',
                charge = -1,
                name = 'Sph%s-CH2-NH2-4H',
                constraints = (
                    FragConstraint(hg = 'Cer', chaintype = 'Sph'),
                ),
                chaintype = 'Sph'
            ), # 267.2330
        # sphingosine long chain base fragments from DH-ceramides
        'Sph_mH_N':
            ChainFragParam(
                plus = 'NH2O2',
                minus = '',
                charge = -1,
                name = 'Sph%s-H',
                constraints = (
                    FragConstraint(hg = 'Cer', chaintype = 'Sph'),
                ),
                chaintype = 'Sph'
            ), # 300.2908
        'Sph_mCH2_mH2O_mH':
            ChainFragParam(
                plus = 'NO',
                minus = 'CH2',
                charge = -1,
                name = 'Sph%s-CH2-H2O-H',
                constraints = (
                    FragConstraint(hg = 'Cer', chaintype = 'Sph'),
                ),
                chaintype = 'Sph'
            ), # 268.2646
        # sphingosine long chain base fragments from d-ceramides
        'Sph_mCH2_mOH':
            ChainFragParam(
                plus = 'NO',
                minus = 'C',
                charge = -1,
                name = 'Sph%s-CH2-OH',
                constraints = (
                    FragConstraint(hg = 'Cer', chaintype = 'Sph'),
                ),
                chaintype = 'Sph'
            ), # 270.2802
        # sphingosine long chain base fragments from sphingomyelin
        'Sph_pPCh_mCH3':
            ChainFragParam(
                plus = 'N2H12O5PC4',
                minus = '',
                charge = -1,
                name = 'Sph%s+PCh-CH3',
                constraints = (
                    FragConstraint(hg = 'SM', chaintype = 'Sph'),
                ),
                chaintype = 'Sph'
            ), # 449.3150
        # sphingosine fragments from sphingosine
        'Sph_mC2H3_mNH2_mOH':
            ChainFragParam(
                plus = 'O',
                minus = 'C2H3',
                charge = -1,
                name = 'Sph%s-C2H3-NH2-OH',
                constraints = (
                    FragConstraint(hg = 'Sph', chaintype = 'Sph'),
                ),
                chaintype = 'Sph'
            ), # 241.2537
        # fatty acyl fragments from ceramide-1-phosphate
        'FA_mH2O':
            ChainFragParam(
                plus = 'O',
                minus = 'H2',
                charge = -1,
                name = 'FA%s-H2O',
                constraints = (
                    FragConstraint(
                        hg = 'Cer',
                        sub = ('1P',),
                        chaintype = 'FA',
                    ),
                ),
                chaintype = 'FA'
            ), # 266.2615
        'FA':
            ChainFragParam(
                plus = 'O2',
                minus = '',
                charge = -1,
                name = 'FA%s-',
                constraints = (
                    FragConstraint(
                        hg = 'Cer',
                        sub = ('1P',),
                        chaintype = 'FA',
                    ),
                ),
                chaintype = 'FA'
            ), # 284.2721
        # sphingosine fragments from sulfohexosyl ceramide
        'Sph_pC6O5H8_pSO3_pH2O':
            ChainFragParam(
                plus = 'C6O11H10S',
                minus = '',
                charge = -1,
                name = 'Sph%s+C6O5H8+SO3+H2O',
                constraints = (
                    FragConstraint(
                        hg = 'Cer',
                        sub = ('SHex',),
                        chaintype = 'Sph',
                    ),
                ),
                chaintype = 'Sph'
            ),
        'Sph_pC6O5H8_pSO3_pCO_pH2O':
            ChainFragParam(
                plus = 'C7O12H10S',
                minus = '',
                charge = -1,
                name = 'Sph%s+C6O5H8+SO3+CO+H2O',
                constraints = (
                    FragConstraint(
                        hg = 'Cer',
                        sub = ('SHex',),
                        chaintype = 'Sph',
                    ),
                ),
                chaintype = 'Sph'
            ),
    }
    
    param_pos = {
        ### neutral losses
        'NLFA':
            ChainFragParam(
                plus = 'O2',
                minus = '',
                charge = 0,
                name = 'NL FA%s',
                constraints = (
                    FragConstraint(chaintype = 'FA'),
                ),
                chaintype = 'FA'
            ),
        # NLFAminusH2O
        'NLFA_mH2O':
            ChainFragParam(
                plus = 'O',
                minus = 'H2',
                charge = 0,
                name = 'NL FA%s-H2O',
                constraints = (
                    FragConstraint(chaintype = 'FA'),
                ),
                chaintype = 'FA'
            ),
        # NLFAplusOH
        'NLFA_pOH':
            ChainFragParam(
                plus = 'O3H',
                minus = '',
                charge = 0,
                name = 'NL FA%s+OH',
                constraints = (
                    FragConstraint(chaintype = 'FA'),
                ),
                chaintype = 'FA'
            ),
        # NLFAplusNH3
        'NLFA_pNH3':
            ChainFragParam(
                plus = 'O2NH3',
                minus = '',
                charge = 0,
                name = 'NL FA%s+NH3',
                constraints = (
                    FragConstraint(family = 'SL', chaintype = 'FA'),
                ),
                chaintype = 'FA'
            ),
        ### positive
        # FAminusO
        'FA_mOH':
            ChainFragParam(
                plus = 'O',
                minus = 'H',
                charge = 1,
                name = 'FA%s-OH',
                constraints = (
                    FragConstraint(chaintype = 'FA'),
                ),
                chaintype = 'FA'
            ), # 267.2682
        # FAplusGlycerol
        'FA_pGlycerol_mOH':
            ChainFragParam(
                plus = 'C3H5O3',
                minus = '',
                charge = 1,
                name = 'FA%s+Glycerol-OH',
                constraints = (
                    FragConstraint(hg = 'PI', chaintype = 'FA'),
                    FragConstraint(hg = 'PG', chaintype = 'FA'),
                    FragConstraint(hg = 'BMP', chaintype = 'FA'),
                    FragConstraint(hg = 'DAG', chaintype = 'FA'),
                    FragConstraint(
                        hg = 'PE',
                        sub = ('Lyso',),
                        chaintype = 'FA',
                    ),
                    FragConstraint(
                        hg = 'PC',
                        sub = ('Lyso',),
                        chaintype = 'FA',
                    ),
                ),
                chaintype = 'FA'
            ),
        # SphingosineBase
        'Sph_pH':
            ChainFragParam(
                plus = 'NH4O2',
                minus = '',
                charge = 1,
                name = 'Sph%s+H',
                constraints = (
                    FragConstraint(
                        hg = 'Sph',
                        chaintype = 'Sph',
                    ),
                    FragConstraint(
                        hg = 'SM',
                        chaintype = 'Sph',
                    ),
                    FragConstraint(
                        hg = 'Cer',
                        chaintype = 'Sph',
                    ),
                    FragConstraint(
                        hg = 'Cer',
                        sub = ('Hex',),
                        chaintype = 'Sph',
                    ),
                ),
                chaintype = 'Sph'
            ), # 302.3053
        # sphingosine long chain base fragments from ceramide-1-phosphate
        'Sph_m2xH2O_pH':
            ChainFragParam(
                plus = 'N',
                minus = '',
                charge = 1,
                name = 'Sph%s-2xH2O+H',
                constraints = (
                    FragConstraint(
                        hg = 'Sph',
                        sph = 'd',
                        chaintype = 'Sph',
                    ),
                    FragConstraint(
                        hg = 'Sph',
                        sph = 'DH',
                        chaintype = 'Sph',
                    ),
                    FragConstraint(
                        hg = 'Sph',
                        sub = ('1P',),
                        chaintype = 'Sph',
                    ),
                    FragConstraint(
                        hg = 'Cer',
                        sub = ('1P',),
                        chaintype = 'Sph',
                    ),
                    FragConstraint(
                        hg = 'SM',
                        chaintype = 'Sph',
                    ),
                    FragConstraint(
                        hg = 'Cer',
                        sph = 'd',
                        chaintype = 'Sph'
                    ),
                    FragConstraint(
                        hg = 'Cer',
                        sph = 'DH',
                        chaintype = 'Sph'
                    ),
                    FragConstraint(
                        hg = 'Cer',
                        sub = ('1P',),
                        chaintype = 'Sph',
                    ),
                    FragConstraint(
                        hg = 'Sph',
                        sub = ('M',),
                        chaintype = 'Sph',
                    ),
                    FragConstraint(
                        hg = 'Sph',
                        sub = ('M2',),
                        chaintype = 'Sph',
                    ),
                ),
                chaintype = 'Sph'
            ), # 266.2842
        # from tCer
        'Sph_m2xH2O_mH':
            ChainFragParam(
                plus = 'N',
                minus = 'H2',
                charge = 1,
                name = 'Sph%s-2xH2O-H',
                constraints = (
                    FragConstraint(
                        hg = 'Cer',
                        sph = 't',
                        chaintype = 'Sph'
                    ),
                    FragConstraint(
                        hg = 'Cer',
                        sub = ('1P',),
                        chaintype = 'Sph',
                    ),
                ),
                chaintype = 'Sph'
            ), # 264.2686
        # fatty acyl fragments from ceramide-1-phosphate
        'FA_pNH_pC2H2_mOH':
            ChainFragParam(
                plus = 'ONC2H2',
                minus = '',
                charge = 1,
                name = 'FA%s+NH+C2H2-OH',
                constraints = (
                    FragConstraint(
                        family = 'SL',
                        chaintype = 'FA'
                    ),
                ),
                chaintype = 'FA'
            ), # 308.2948
        'FA_pNH2_mO':
            ChainFragParam(
                plus = 'ONH2',
                minus = '',
                charge = 1,
                name = 'FA%s+NH2-O',
                constraints = (
                    FragConstraint(
                        hg = 'Cer',
                        sph = 'd',
                        chaintype = 'FA',
                    ),
                    FragConstraint(
                        hg = 'Cer',
                        sph = 'DH',
                        chaintype = 'FA',
                    ),
                    FragConstraint(
                        hg = 'Cer',
                        sub = ('1P',),
                        chaintype = 'FA',
                    ),
                    FragConstraint(
                        hg = 'Cer',
                        sub = ('PE',),
                        chaintype = 'FA',
                    ),
                ),
                chaintype = 'FA'
            ), # 284.2948
        # sphingosine long chain base fragments from hexosyl-t-ceramide
        # also from tCer
        'Sph_mH_P':
            ChainFragParam(
                plus = 'NH2O2',
                minus = '',
                charge = 1,
                name = 'Sph%s-H',
                constraints = (
                    FragConstraint(
                        hg = 'Sph',
                        sph = 't',
                        chaintype = 'Sph',
                    ),
                    FragConstraint(
                        hg = 'Cer',
                        sph = 't',
                        chaintype = 'Sph',
                    ),
                    FragConstraint(
                        hg = 'Cer',
                        sub = ('Hex',),
                        chaintype = 'Sph',
                    ),
                ),
                chaintype = 'Sph'
            ), # 300.2897
        # also from tCer & tSph
        'Sph_pH2O_mH':
            ChainFragParam(
                plus = 'NH4O3',
                minus = '',
                charge = 1,
                name = 'Sph%s+H2O-H',
                constraints = (
                    FragConstraint(
                        hg = 'Sph',
                        sph = 't',
                        chaintype = 'Sph',
                    ),
                    FragConstraint(
                        hg = 'Cer',
                        sub = ('Hex',),
                        chaintype = 'Sph',
                    ),
                ),
                chaintype = 'Sph'
            ), # 318.3003
        # also from tCer
        'Sph_mH2O_mH':
            ChainFragParam(
                plus = 'NO',
                minus = '',
                charge = 1,
                name = 'Sph%s-H2O-H',
                constraints = (
                    FragConstraint(
                        hg = 'Cer',
                        sph = 't',
                        chaintype = 'Sph',
                    ),
                    FragConstraint(
                        hg = 'Cer',
                        sub = ('Hex',),
                        chaintype = 'Sph',
                    ),
                    FragConstraint(
                        hg = 'Sph',
                        sph = 't',
                        chaintype = 'Sph',
                    ),
                ),
                chaintype = 'Sph'
            ), # 282.2791
        'Sph_mH2O_pH':
            ChainFragParam(
                plus = 'NOH2',
                minus = '',
                charge = 1,
                name = 'Sph%s-H2O+H',
                constraints = (
                    FragConstraint(
                        hg = 'Cer',
                        sph = 'DH',
                        chaintype = 'Sph',
                    ),
                    FragConstraint(
                        hg = 'Sph',
                        sub = ('1P',),
                        chaintype = 'Sph',
                    ),
                    FragConstraint(
                        hg = 'Sph',
                        sub = ('M',),
                        chaintype = 'Sph',
                    ),
                    FragConstraint(
                        hg = 'Sph',
                        sub = ('M2',),
                        chaintype = 'Sph',
                    ),
                    FragConstraint(
                        hg = 'Sph',
                        sph = 'd',
                        chaintype = 'Sph',
                    ),
                    FragConstraint(
                        hg = 'Sph',
                        sph = 'DH',
                        chaintype = 'Sph',
                    ),
                    FragConstraint(
                        hg = 'Cer',
                        sph = 'd',
                        chaintype = 'Sph',
                    ),
                ),
                chaintype = 'Sph'
            ), # 284.2948
        # sphingosine fragments from N-(n)methyl-safingols
        'Sph_mH2O_pCH3':
            ChainFragParam(
                plus = 'NOCH4',
                minus = '',
                charge = 1,
                name = 'Sph%s-H2O+CH3',
                constraints = (
                    FragConstraint(
                        hg = 'Sph',
                        sub = ('M',),
                        chaintype = 'Sph',
                    ),
                    FragConstraint(
                        hg = 'Sph',
                        sub = ('M2',),
                        chaintype = 'Sph',
                    ),
                ),
                chaintype = 'Sph'
            ), # 298.3104
        'Sph_mH2O_p2xCH3_mH':
            ChainFragParam(
                plus = 'NOC2H6',
                minus = '',
                charge = 1,
                name = 'Sph%s-H2O+2xCH3-H',
                constraints = (
                    FragConstraint(
                        hg = 'Sph',
                        sub = ('M2',),
                        chaintype = 'Sph',
                    ),
                ),
                chaintype = 'Sph'
            ), # 312.3261
        # sphinosine fragments from keto-sphingosine & tSph
        'Sph_mC_mH2O_mH':
            ChainFragParam(
                plus = 'NO',
                minus = 'C',
                charge = 1,
                name = 'Sph%s-C-H2O-H',
                constraints = (
                    FragConstraint(
                        hg = 'Sph',
                        sph = 'k',
                        chaintype = 'Sph',
                    ),
                    FragConstraint(
                        hg = 'Sph',
                        sph = 't',
                        chaintype = 'Sph',
                    ),
                ),
                chaintype = 'Sph'
            ), # 270.2791
        'Sph_mC_mO_mH2O_mH':
            ChainFragParam(
                plus = 'N',
                minus = 'C',
                charge = 1,
                name = 'Sph%s-C-O-H2O-H',
                constraints = (
                    FragConstraint(
                        hg = 'Cer',
                        sph = 'd',
                        chaintype = 'Sph',
                    ),
                    FragConstraint(
                        hg = 'Cer',
                        sph = 'DH',
                        chaintype = 'Sph',
                    ),
                    FragConstraint(
                        hg = 'Sph',
                        sph = 'k',
                        chaintype = 'Sph',
                    ),
                    FragConstraint(
                        hg = 'Sph',
                        sph = 'd',
                        chaintype = 'Sph',
                    ),
                    FragConstraint(
                        hg = 'Sph',
                        sph = 'DH',
                        chaintype = 'Sph',
                    ),
                    FragConstraint(
                        hg = 'Sph',
                        sph = 't',
                        chaintype = 'Sph',
                    ),
                ),
                chaintype = 'Sph'
            ), # 254.2842
        'Sph_mC_mO_mH2O_mNH':
            ChainFragParam(
                plus = '',
                minus = 'C',
                charge = 1,
                name = 'Sph%s-C-O-H2O-NH',
                constraints = (
                    FragConstraint(
                        hg = 'Cer',
                        sph = 'DH',
                        chaintype = 'Sph',
                    ),
                    FragConstraint(
                        hg = 'Sph',
                        sph = 'DH',
                        chaintype = 'Sph',
                    ),
                    FragConstraint(
                        hg = 'Sph',
                        sph = 't',
                        chaintype = 'Sph',
                    ),
                ),
                chaintype = 'Sph'
            ), # 240.2812
        'Sph_mC2H2_mO_mH2O_mNH':
            ChainFragParam(
                plus = '',
                minus = 'C2H2',
                charge = 1,
                name = 'Sph%s-C2H2-O-H2O-NH',
                constraints = (
                    FragConstraint(
                        hg = 'Sph',
                        sph = 't',
                        chaintype = 'Sph',
                    ),
                ),
                chaintype = 'Sph'
            ), # 226.2655
        'Sph_mNH2_mH2O_m2H':
            ChainFragParam(
                plus = 'O',
                minus = 'H3',
                charge = 1,
                name = 'Sph%s-C-H2O-2H',
                constraints = (
                    FragConstraint(
                        hg = 'Sph',
                        sph = 'k',
                        chaintype = 'Sph',
                    ),
                ),
                chaintype = 'Sph'
            ), # 265.2526
        'Sph_mC_m2xH2O':
            ChainFragParam(
                plus = 'N',
                minus = 'CH2',
                charge = 1,
                name = 'Sph%s-C-2xH2O',
                constraints = (
                    FragConstraint(
                        hg = 'Sph',
                        sph = 'k',
                        chaintype = 'Sph',
                    ),
                    FragConstraint(
                        hg = 'Cer',
                        sph = 't',
                        chaintype = 'Sph',
                    ),
                    FragConstraint(
                        hg = 'Sph',
                        sph = 't',
                        chaintype = 'Sph',
                    ),
                ),
                chaintype = 'Sph'
            ), # 252.2686
        # sphingosine fragments from d-ceramide
        'Sph_mO_pH':
            ChainFragParam(
                plus = 'NOH4',
                minus = '',
                charge = 1,
                name = 'Sph%s-O+H',
                constraints = (
                    FragConstraint(
                        hg = 'Sph',
                        sph = 'd',
                        chaintype = 'Sph',
                    ),
                ),
                chaintype = 'Sph'
            ), # 286.3104
    }
    
    def __init__(self):
        
        mod = sys.modules[__name__]
        
        for ionmode in ('neg', 'pos'):
            
            param = getattr(self, 'param_%s' % ionmode)
            
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
                        '        ionmode = \'%s\',\n'
                        '        attrs = {\n'
                        '            \'constraints\': %s,\n'
                        '            \'chaintype\': \'%s\'\n'
                        '        },\n'
                        '        **kwargs\n'
                        '    )\n'
                    ) % (
                        par.plus,
                        par.minus,
                        par.name,
                        par.charge,
                        ionmode,
                        str(par.constraints),
                        par.chaintype
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
                    ('chaintype', par.chaintype),
                    ('constraints', par.constraints),
                    ('name', par.name % ''),
                    ('ionmode', ionmode)
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
