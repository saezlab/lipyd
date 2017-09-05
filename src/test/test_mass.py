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

import unittest
import emese.mass as mass

class TestFormula(unittest.TestCase):
    
    def test_calc_weight(self):
        
        ethanol = mass.Formula('C2H5OH')
        self.assertTrue(abs(ethanol.weight - 46.04186481376) < 0.0000001)
    
    def test_bind(self):
        
        ethanol = mass.Formula('C2H5OH')
        aceticacid = mass.Formula('CH3COOH')
        diethylester = ethanol.bind(aceticacid)
        self.assertEqual(diethylester.formula, 'C4H8O2')
    
    def test_divide(self):
        
        ethanol = mass.Formula('C2H5OH')
        aceticacid = mass.Formula('CH3COOH')
        diethylester = ethanol.bind(aceticacid)
        ethanol2, aceticacid2 = diethylester.divide(ethanol)
        
        self.assertEqual(ethanol.formula, ethanol2.formula)
        self.assertEqual(aceticacid.formula, aceticacid2.formula)
    
    def test_add(self):
        
        aceticacid  = mass.Formula('CH3COOH')
        acetate     = aceticacid - mass.Formula('H', charge = 1)
        acetate += mass.Formula('H', charge = 1)
        self.assertTrue(abs(acetate.weight - aceticacid.weight) < 0.0000001)
