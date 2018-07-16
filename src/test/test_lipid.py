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
        903.6717059967399
    )
]


class TestLipid(object):
    
    @pytest.mark.parametrize('clsname, args, mass', gpl)
    def test_glycerophospholipids(self, clsname, args, mass):
        
        cls = getattr(lipyd.lipid, clsname)
        lip = list(cls(**args))[0]
        
        assert abs(lip.mass - mass) < 0.000001
