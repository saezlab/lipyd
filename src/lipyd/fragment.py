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

from __future__ import print_function
from future.utils import iteritems
from past.builtins import xrange, range, reduce

import re
import imp

import lipyd.mass as mass

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


class FattyFragment(mass.MassBase, AdductCalculator):
    
    def __init__(self, charge, c = 3, unsat = 0,
        minus = [], plus = [], isotope = 0, name = None, hg = []):
        self.c = c
        self.unsat = unsat
        self.minus = minus
        self.plus = plus
        self.hg = hg
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


class LysoPEAlkenyl(FattyFragment):
    """
    from massbank.jp:
    [lyso PE(alkenyl-18:0,-)]- 464.3140997565 -417 C23H47NO6P-
    """
    
    def __init__(self, c, unsat = 0,
        minus = [], plus = [], isotope = 0,
        charge = -1):
        self.counts = {
            'C': 5,
            'H': 11,
            'O': 4,
            'N': 1,
            'P': 1
        }
        super(LysoPEAlkenyl, self).__init__(
            charge = charge,
            c = c,
            unsat = unsat,
            minus = minus,
            plus = plus,
            isotope = isotope,
            name = 'Lyso-PE-alkenyl',
            hg = ['PE']
        )


class LysoPE(FattyFragment):
    """
    from massbank.jp:
    [lyso PE(18:0,-)]- 480.3090143786 -476 C23H47NO7P-
    """
    
    def __init__(self, c, unsat = 0,
        minus = [], plus = [], isotope = 0,
        charge = -1):
        self.counts = {
            'C': 5,
            'H': 11,
            'O': 5,
            'N': 1,
            'P': 1
        }
        super(LysoPE, self).__init__(
            charge = charge,
            c = c,
            unsat = unsat,
            minus = minus,
            plus = plus,
            isotope = isotope,
            name = 'Lyso-PE',
            hg = ['PE']
        )


class LysoPEAlkyl(FattyFragment):
    
    def __init__(self, c, unsat = 0,
        minus = [], plus = [], isotope = 0,
        charge = -1):
        self.counts = {
            'C': 5,
            'H': 13,
            'O': 4,
            'N': 1,
            'P': 1
        }
        super(LysoPEAlkyl, self).__init__(
            charge = charge,
            c = c,
            unsat = unsat,
            minus = minus,
            plus = plus,
            isotope = isotope,
            name = 'Lyso-PE-alkyl',
            hg = ['PE']
        )

class LysoPCAlkenyl(FattyFragment):
    """
    from massbank.jp:
    [lyso PC(alkenyl-18:0,-)]- 492.3453998849 -436 C25H51NO6P-
    """
    
    def __init__(self, c, unsat = 0,
        minus = [], plus = [], isotope = 0,
        charge = -1):
        self.counts = {
            'C': 7,
            'H': 15,
            'O': 4,
            'N': 1,
            'P': 1
        }
        super(LysoPCAlkenyl, self).__init__(
            charge = charge,
            c = c,
            unsat = unsat,
            minus = minus,
            plus = plus,
            isotope = isotope,
            name = 'Lyso-PC-alkenyl',
            hg = ['PC']
        )


class LysoPCAlkyl(FattyFragment):
    """
    from massbank.jp:
    [lyso PC(alkyl-18:0,-)]- 494.3610499491 -143 C25H53NO6P-
    """
    
    def __init__(self, c, unsat = 0,
        minus = [], plus = [], isotope = 0,
        charge = -1):
        self.counts = {
            'C': 7,
            'H': 17,
            'O': 4,
            'N': 1,
            'P': 1
        }
        super(LysoPCAlkyl, self).__init__(
            charge = charge,
            c = c,
            unsat = unsat,
            minus = minus,
            plus = plus,
            isotope = isotope,
            name = 'Lyso-PC-alkyl',
            hg = ['PC']
        )


class LysoPC(FattyFragment):
    
    # from massbank.jp:
    # [lyso PC(18:0,-)]- 508.340314507 -373 C25H51NO7P-
    
    def __init__(self, c, unsat = 0,
        minus = [], plus = [], isotope = 0,
        charge = -1):
        self.counts = {
            'C': 7,
            'H': 15,
            'O': 5,
            'N': 1,
            'P': 1
        }
        super(LysoPC, self).__init__(
            charge = charge,
            c = c,
            unsat = unsat,
            minus = minus,
            plus = plus,
            isotope = isotope,
            name = 'Lyso-PC',
            hg = ['PC']
        )


class LysoPS(FattyFragment):
    """
    from massbank.jp:
    [lyso PS(18:0,-)]- 437.2668152129 -358 C21H42O7P-
    """
    
    def __init__(self, c, unsat = 0,
        minus = [], plus = [], isotope = 0,
        charge = -1):
        self.counts = {
            'C': 3,
            'H': 6,
            'O': 5,
            'P': 1
        }
        super(LysoPS, self).__init__(
            charge = charge,
            c = c,
            unsat = unsat,
            minus = minus,
            plus = plus,
            isotope = isotope,
            name = 'Lyso-PS',
            hg = ['PS']
        )


class LysoPSAlkenyl(FattyFragment):
    
    def __init__(self, c, unsat = 0,
        minus = [], plus = [], isotope = 0,
        charge = -1):
        self.counts = {
            'C': 3,
            'H': 6,
            'O': 4,
            'P': 1
        }
        super(LysoPSAlkenyl, self).__init__(
            charge = charge,
            c = c,
            unsat = unsat,
            minus = minus,
            plus = plus,
            isotope = isotope,
            name = 'Lyso-PS-alkenyl',
            hg = ['PS']
        )


class LysoPSAlkyl(FattyFragment):
    
    def __init__(self, c, unsat = 0,
        minus = [], plus = [], isotope = 0,
        charge = -1):
        self.counts = {
            'C': 3,
            'H': 8,
            'O': 4,
            'P': 1
        }
        super(LysoPSAlkyl, self).__init__(
            charge = charge,
            c = c,
            unsat = unsat,
            minus = minus,
            plus = plus,
            isotope = isotope,
            name = 'Lyso-PS-alkyl',
            hg = ['PS']
        )


class LysoPI(FattyFragment):
    """
    from massbank:
    [lyso PI(-,18:0)]- 599.3196386444 -432 C27H52O12P-
    """
    
    def __init__(self, c, unsat = 0,
        minus = [], plus = [], isotope = 0,
        charge = -1):
        self.counts = {
            'C': 9,
            'H': 16,
            'O': 10,
            'P': 1
        }
        super(LysoPI, self).__init__(
            charge = charge,
            c = c,
            unsat = unsat,
            minus = minus,
            plus = plus,
            isotope = isotope,
            name = 'Lyso-PI',
            hg = ['PI']
        )


class LysoPIAlkyl(FattyFragment):
    """
    from massbank.jp:
    [lyso PI(alkyl-16:1,-)-H2O]- 537.2828592076 -228 C25H46O10P-
    from this, derived 18:0-:
    [lyso PI(alkyl-18:0,-)]- 585.3403740858 -228 C27H54O11P-
    """
    
    def __init__(self, c, unsat = 0,
        minus = [], plus = [], isotope = 0,
        charge = -1):
        self.counts = {
            'C': 9,
            'H': 18,
            'O': 9,
            'P': 1
        }
        super(LysoPIAlkyl, self).__init__(
            charge = charge,
            c = c,
            unsat = unsat,
            minus = minus,
            plus = plus,
            isotope = isotope,
            name = 'Lyso-PI-alkyl',
            hg = ['PI']
        )


class LysoPG(FattyFragment):
    """
    from massbank:
    [lyso PG(18:0,-)]- 511.3035946497 -495 C24H48O9P-
    """
    
    def __init__(self, c, unsat = 0,
        minus = [], plus = [], isotope = 0,
        charge = -1):
        self.counts = {
            'C': 6,
            'H': 12,
            'O': 7,
            'P': 1
        }
        super(LysoPG, self).__init__(
            charge = charge,
            c = c,
            unsat = unsat,
            minus = minus,
            plus = plus,
            isotope = isotope,
            name = 'Lyso-PG',
            hg = ['PG']
        )


class LysoPGAlkyl(FattyFragment):
    """
    from massbank:
    [lyso PG(18:0,-)]- 511.3035946497 -495 C24H48O9P-
    """
    
    def __init__(self, c, unsat = 0,
        minus = [], plus = [], isotope = 0,
        charge = -1):
        self.counts = {
            'C': 6,
            'H': 14,
            'O': 6,
            'P': 1
        }
        super(LysoPGAlkyl, self).__init__(
            charge = charge,
            c = c,
            unsat = unsat,
            minus = minus,
            plus = plus,
            isotope = isotope,
            name = 'Lyso-PG-alkyl',
            hg = ['PG']
        )


class LysoPA(FattyFragment):
    """
    from Characterization of Phospholipid 
    Molecular Species by Means of HPLC-Tandem Mass Spectrometry:
    [lyso PA(18:1,-)]- 435.2 C21H40O7P-
    """
    
    def __init__(self, c, unsat = 0,
        minus = [], plus = [], isotope = 0,
        charge = -1):
        self.counts = {
            'C': 3,
            'H': 6,
            'O': 5,
            'P': 1
        }
        super(LysoPA, self).__init__(
            charge = charge,
            c = c,
            unsat = unsat,
            minus = minus,
            plus = plus,
            isotope = isotope,
            name = 'Lyso-PA',
            hg = ['PA']
        )


class LysoPAAlkyl(FattyFragment):
    """
    from Characterization of Phospholipid 
    Molecular Species by Means of HPLC-Tandem Mass Spectrometry:
    [lyso PA(18:1,-)]- 435.2 C21H40O7P-
    """
    
    def __init__(self, c, unsat = 0,
        minus = [], plus = [], isotope = 0,
        charge = -1):
        self.counts = {
            'C': 3,
            'H': 8,
            'O': 4,
            'P': 1
        }
        super(LysoPAAlkyl, self).__init__(
            charge = charge,
            c = c,
            unsat = unsat,
            minus = minus,
            plus = plus,
            isotope = isotope,
            name = 'Lyso-PA-alkyl',
            hg = ['PA']
        )

class CerFA(FattyFragment):
    """
    https://metlin.scripps.edu/metabo_info.php?molid=6214
    [Cer-FA(C8:0)]- 168.1382904 C8H14O1N1-
    """
    
    def __init__(self, c, unsat = 0,
        minus = [], plus = [], isotope = 0,
        charge = -1):
        self.counts = {
            'C': 2,
            'H': 2,
            'O': -1,
            'N': 1
        }
        super(CerFA, self).__init__(
            charge = charge,
            c = c,
            unsat = unsat,
            minus = minus,
            plus = plus,
            isotope = isotope,
            name = 'CerFA',
            hg = ['Cer']
        )


class CerFAminusC2H5N(FattyFragment):
    
    # 209 at C14:0 FA, 263 at C18:1
    
    def __init__(self, c, unsat = 0,
        minus = [], plus = [], isotope = 0,
        charge = -1):
        self.counts = {
            'C': 0,
            'H': -3,
            'O': -1
        }
        super(CerFAminusC2H5N, self).__init__(
            charge = charge,
            c = c,
            unsat = unsat,
            minus = minus,
            plus = plus,
            isotope = isotope,
            name = 'CerFA-C2N',
            hg = ['Cer']
        )


class CerFAminusC(FattyFragment):
    
    # 226 at C14:0 FA, 280 at C18:1
    
    def __init__(self, c, unsat = 0,
        minus = [], plus = [], isotope = 0,
        charge = -1):
        self.counts = {
            'C': 0,
            'H': 0,
            'O': -1,
            'N': 1
        }
        super(CerFAminusC, self).__init__(
            charge = charge,
            c = c,
            unsat = unsat,
            minus = minus,
            plus = plus,
            isotope = isotope,
            name = 'CerFA-C',
            hg = ['Cer']
        )


class CerFAminusN(FattyFragment):
    
    # 227 at C14:0 FA, 281 at C18:1
    
    def __init__(self, c, unsat = 0,
        minus = [], plus = [], isotope = 0,
        charge = -1):
        self.counts = {
            'C': 0,
            'H': -1,
            'O': 0,
        }
        super(CerFAminusN, self).__init__(
            charge = charge,
            c = c,
            unsat = unsat,
            minus = minus,
            plus = plus,
            isotope = isotope,
            name = 'CerFA-N',
            hg = ['Cer']
        )


class CerSphiMinusN(FattyFragment):
    
    def __init__(self, c, unsat = 0,
        minus = [], plus = [], isotope = 0,
        charge = -1):
        self.counts = {
            'C': -2,
            'H': -5,
            'O': -1
        }
        super(CerSphiMinusN, self).__init__(
            charge = charge,
            c = c,
            unsat = unsat,
            minus = minus,
            plus = plus,
            isotope = isotope,
            name = 'CerSphi-N',
            hg = ['Cer']
        )


class CerSphiMinusNO(FattyFragment):
    
    def __init__(self, c, unsat = 0,
        minus = [], plus = [], isotope = 0,
        charge = -1):
        self.counts = {
            'C': 0,
            'H': -3,
            'O': -1
        }
        super(CerSphiMinusNO, self).__init__(
            charge = charge,
            c = c,
            unsat = unsat,
            minus = minus,
            plus = plus,
            isotope = isotope,
            name = 'CerSphi-N-H2O',
            hg = ['Cer']
        )


class CerSphi(FattyFragment):
    
    def __init__(self, c, unsat = 0,
        minus = [], plus = [], isotope = 0,
        charge = -1):
        self.counts = {
            'C': -2,
            'H': -4,
            'O': 0,
            'N': 1
        }
        super(CerSphi, self).__init__(
            charge = charge,
            c = c,
            unsat = unsat,
            minus = minus,
            plus = plus,
            isotope = isotope,
            name = 'CerSphi',
            hg = ['Cer']
        )


class FAminusH(FattyFragment):
    
    def __init__(self, c, unsat = 0, isotope = 0, **kwargs):
        super(FAminusH, self).__init__(charge = -1, c = c, unsat = unsat,
            minus = ['H'], isotope = isotope, name = 'FA')


class FAAlkylminusH(FattyFragment):
    """
    18:0 = 267.2693393
    """
    
    def __init__(self, c, unsat = 1, isotope = 0, **kwargs):
        self.counts = {
            'O': -1,
            'H': 2
        }
        super(FAAlkylminusH, self).__init__(charge = -1, c = c, unsat = unsat,
            minus = ['H'], isotope = isotope, name = 'FA-alkyl')


class NLFA(FattyFragment):
    
    def __init__(self, c, unsat = 0, isotope = 0, **kwargs):
        super(NLFA, self).__init__(charge = 0, c = c, unsat = unsat,
            isotope = isotope, name = 'NL FA')


class NLFAminusH2O(FattyFragment):
    
    def __init__(self, c, unsat = 0, isotope = 0, **kwargs):
        super(NLFAminusH2O, self).__init__(charge = 0, c = c, unsat = unsat,
            minus = ['H2O'], isotope = isotope, name = 'NL FA')


class NLFAplusOH(FattyFragment):
    
    def __init__(self, c, unsat = 0, isotope = 0, **kwargs):
        super(NLFAplusOH, self).__init__(charge = 0, c = c, unsat = unsat,
            plus = ['OH'], isotope = isotope, name = 'NL FA')


class NLFAplusNH3(FattyFragment):
    
    def __init__(self, c, unsat = 0, isotope = 0, **kwargs):
        super(NLFAplusNH3, self).__init__(charge = 0, c = c, unsat = unsat,
            plus = ['NH3'], isotope = isotope, name = 'NL FA')


class FAminusO(FattyFragment):
    
    def __init__(self, c, unsat = 0, isotope = 0, **kwargs):
        super(FAminusO, self).__init__(charge = 1, c = c, unsat = unsat,
            minus = ['O', 'H'], isotope = isotope, name = 'FA')


class FAplusGlycerol(FattyFragment):
    
    def __init__(self, c, unsat = 0, isotope = 0,
        minus = [], plus = [], **kwargs):
        self.counts = {
            'C': 3,
            'H': 5,
            'O': 1
        }
        super(FAplusGlycerol, self).__init__(charge = 1, c = c, unsat = unsat,
            minus = minus, plus = plus,
            isotope = isotope,
            name = 'FA+G',
            hg = ['PG', 'BMP', 'DAG', 'LysoPE', 'LysoPC']
        )


class SphingosineBase(FattyFragment):
    
    def __init__(self, c, unsat = 0, isotope = 0,
        minus = [], plus = [], **kwargs):
        self.counts = {
            'N': 1,
            'H': 4
        }
        super(SphingosineBase, self).__init__(charge = 1, c = c, unsat = unsat,
            isotope = isotope, name = 'Sphingosine',
            minus = minus, plus = plus,
            hg = ['Sph', 'SM', 'Cer', 'HexCer'])


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
