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

import lipyd.formula as formula


class TestFormula(object):
    
    def test_calc_mass(self):
        
        ethanol = formula.Formula('C2H5OH')
        
        assert abs(ethanol.mass - 46.04186481376) < 0.0000001
    
    def test_bind(self):
        
        ethanol = formula.Formula('C2H5OH')
        aceticacid = formula.Formula('CH3COOH')
        diethylester = ethanol.bind(aceticacid)
        
        assert diethylester.formula == 'C4H8O2'
    
    def test_bind_split(self):
        
        ethanol = formula.Formula('C2H5OH')
        aceticacid = formula.Formula('CH3COOH')
        diethylester = ethanol.bind(aceticacid)
        ethanol2, aceticacid2 = diethylester.split(ethanol)
        
        assert ethanol.formula == ethanol2.formula
        
        assert aceticacid.formula == aceticacid2.formula
    
    def test_add(self):
        
        aceticacid  = formula.Formula('CH3COOH')
        acetate     = aceticacid - formula.Formula('H', charge = 1)
        
        assert acetate.charge == -1
        
        acetate    += formula.Formula('H', charge = 1)
        
        assert abs(acetate.mass - aceticacid.mass) < 0.0000001
