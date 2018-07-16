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
    (
        'PI',
        {'fa_args': {'c': 38, 'u': 5}, 'sn2_fa_args': {'c': 18, 'u': 3}},
        1130.7762306419602,
        'PI(38:5/18:3)'
    ),
    (
        'PGP',
        {'fa_args': {'c': 16, 'u': 2}, 'sn2_fa_args': {'c': 26, 'u': 5}},
        928.5230667049199,
        'PGP(16:2/26:5)'
    ),
    (
        'PC',
        {'fa_args': {'c': 18, 'u': 3}, 'sn2_fa_args': {'c': 34, 'u': 4}},
        999.76560638386,
        'PC(18:3/34:4)'
    ),
    (
        'EtherPA',
        {'fa_args': {'c': 15, 'u': 0}, 'sn2_fa_args': {'c': 30, 'u': 4}},
        808.6345922110401,
        'PA(O-15:0/30:4)'
    ),
    (
        'PS',
        {'fa_args': {'c': 34, 'u': 6}, 'sn2_fa_args': {'c': 28, 'u': 5}},
        1133.80238581782,
        'PS(34:6/28:5)'
    ),
    (
        'PG',
        {'fa_args': {'c': 20, 'u': 4}, 'sn2_fa_args': {'c': 34, 'u': 5}},
        1012.71323645876,
        'PG(20:4/34:5)'
    ),
    (
        'PIP',
        {'fa_args': {'c': 32, 'u': 5}, 'sn2_fa_args': {'c': 13, 'u': 0}},
        1062.61736101716,
        'PIP(32:5/13:0)'
    ),
    (
        'PIP2',
        {'fa_args': {'c': 16, 'u': 1}, 'sn2_fa_args': {'c': 38, 'u': 5}},
        1266.70889242468,
        'PIP2(16:1/38:5)'
    ),
    (
        'PIP3',
        {'fa_args': {'c': 18, 'u': 1}, 'sn2_fa_args': {'c': 36, 'u': 5}},
        1346.6752233160398,
        'PIP3(18:1/36:5)'
    )
]

gl = [
    (
        'MonoalkylGlycerol',
        {'fa_args': {'c': 18, 'u': 2}},
        340.297745151,
        'MAG(O-18:2)'
    ),
    (
        'MonoalkylMonoacylGlycerol',
        {'fa_args': {'c': 20, 'u': 1}, 'sn2_fa_args': {'c': 24, 'u': 6}},
        708.6056610616,
        'DAG(O-20:1/24:6)'
    ),
    (
        'MonoalkylDiacylGlycerol',
        {
            'fa_args': {'c': 18, 'u': 0},
            'sn2_fa_args': {'c': 17, 'u': 0},
            'sn3_fa_args': {'c': 22, 'u': 0}
        },
        918.89792690768,
        'TAG(O-18:0/17:0/22:0)'
    ),
    (
        'TriacylGlycerol',
        {
            'fa_args': {'c': 18, 'u': 1},
            'sn2_fa_args': {'c': 26, 'u': 0},
            'sn3_fa_args': {'c': 16, 'u': 1}
        },
        970.89284152788,
        'TAG(18:1/26:0/16:1)'
    ),
    (
        'DiacylGlycerol',
        {
            'fa_args': {'c': 18, 'u': 4},
            'sn2_fa_args': {'c': 24, 'u': 5},
        },
        690.5223253592001,
        'DAG(18:4/24:5)'
    )
]


class TestLipid(object):
    
    @pytest.mark.parametrize('clsname, args, mass, name', gpl)
    def test_glycerophospholipids(self, clsname, args, mass, name):
        
        cls = getattr(lipyd.lipid, clsname)
        lip = list(cls(**args))[0]
        
        assert abs(lip.mass - mass) < 0.000001
        assert lip.name == name
    
    @pytest.mark.parametrize('clsname, args, mass, name', gl)
    def test_glycerolipids(self, clsname, args, mass, name):
        
        cls = getattr(lipyd.lipid, clsname)
        lip = list(cls(**args))[0]
        
        assert abs(lip.mass - mass) < 0.000001
        assert lip.name == name
