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

import pytest

import lipyd.fragment


standards = {
    'FA_mH': 283.26425394940924,
    'LysoPE': 480.3095634402292,
    'LysoPC': 508.34086356926923,
    'LysoPA': 437.26736427472923,
    'LysoPG': 511.30414370868925,
    'LysoPI': 599.3201876983293,
    'FAL_mH': 269.28498939372923,
    'LysoPA_mH2O': 419.2567995900092,
    'FA_mH2O_mH': 265.2536892646893,
    'FA_mO_pNH2': 282.2802383656693,
    'FA_mO_pC2H2NH2': 308.2958884301893,
    'Sph_mC2H4_mNH2_mH2O': 239.23803920016925,
    'Sph_mNH2_mH2O_m2H': 265.25259210487076,
    'Sph_mC2H4_m3H': 270.2438528568293,
    'NLFA': 284.27153040176,
    'NLFA_mH2O': 266.26096571704,
    'NLFA_pOH': 301.27427005422,
    'NLFA_pNH3': 301.29807950274,
    'FA_mOH': 267.2682421693907,
    'FA_pGlycerol_mOH': 341.3050216033507,
    'Sph_pH': 302.30535595509076
}


class TestFragment(object):
    
    def test_fatty_fragment(self):
        
        ff = list(lipyd.fragment.FattyFragment(
            c = 18, u = 0, charge = -1, head = 'OO', minus = 'H'
        ))[0]
        
        assert ff.formula == 'C18H35O2'
        assert ff.charge  == -1
        assert abs(ff.mass - 283.26425394940924) < 0.0000001
    
    @pytest.mark.parametrize('name, refmass', standards.items())
    def test_fragment_classes(self, name, refmass):
        
        toalkyl = (
            lipyd.fragment.mass.MassBase('O').mass -
            lipyd.fragment.mass.MassBase('H2').mass
        )
        
        for name in lipyd.fragment.fattyfragments:
            
            if name in standards:
                
                assert hasattr(lipyd.fragment, name)
                
                cls = getattr(lipyd.fragment, name)
                ff = list(cls(c = 18, u = 0))[0]
                
                assert abs(ff.mass - standards[name]) < 0.0000001
                
                if name.startswith('Lyso') and not name.endswith('Alkyl'):
                    
                    aname = '%sAlkyl' % name
                    
                    if not hasattr(lipyd.fragment, aname):
                        continue
                    
                    acls = getattr(lipyd.fragment, aname)
                    aff = list(acls(c = 18, u = 0))[0]
                    
                    assert abs(aff.mass + toalkyl - ff.mass) < 0.0000001
    
    def test_fragment_names(self):
        
        fr = list(lipyd.fragment.Sph_mNH2_mH2O_m2H(c = 18, u = 1))[0]
        
        assert fr.name == '[Sph(18:1)-NH2-H2O-2H]+'
    
    def test_fatty_fraglines(self):
        
        sample_line = [
            '[Sph(18:1)-NH2-H2O-2H]+', '', '',
            18, 1, -1, 'Sph-NH2-H2O-2H'
        ]
        
        fl = list(
            lipyd.fragment.Sph_mNH2_mH2O_m2H(
                c = (18), u = (1)
            ).iterfraglines()
        )
        
        for real, sample in zip(fl[1:], sample_line):
            
            assert real == sample
