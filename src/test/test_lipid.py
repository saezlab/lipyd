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

import pytest

import lipyd.lipid

gpl = [
    (
        'PE',
        {'fa_args': {'c': 26, 'u': 0}, 'sn2_fa_args': {'c': 22, 'u': 6}},
        903.6717059967399,
        'PE(26:0/22:6)'
    ),
    (
        'EtherPE',
        {'fa_args': {'c': 16, 'u': 0}, 'sn2_fa_args': {'c': 28, 'u': 0}},
        845.7237415701001,
        'PE(O-16:0/28:0)'
    ),
    (
        'LysoPE',
        {'fa_args': {'c': 18, 'u': 1}},
        479.30118982806,
        'LysoPE(18:1)'
    ),
]


class TestLipid(object):
    
    @pytest.mark.parametrize('clsname, args, mass, name', gpl)
    def test_glycerophospholipids(self, clsname, args, mass, name):
        
        cls = getattr(lipyd.lipid, clsname)
        lip = list(cls(**args))[0]
        
        assert abs(lip.mass - mass) < 0.000001
        assert lip.name == name
