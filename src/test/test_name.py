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

import lipyd.name

swl_names = {
    'Phosphatidylcholine(36:1)':
        ('PC', ['', 36, 1], [], None),
    'Phosphatidylcholine(18:1/18:0)':
        ('PC', ['', 36, 1], [('', 18, 1), ('', 18, 0)], None)
}

lmp_names = {
    'Cer(t18:1/18:1)':
        ('tCer', ['t', 36, 2], [('t', 18, 1), ('', 18, 1)], None),
    'Cer(d18:1/18:1)':
        ('dCer', ['d', 36, 2], [('d', 18, 1), ('', 18, 1)], None),
    'Cer(d18:0/18:1)':
        ('DHCer', ['d', 36, 1], [('d', 18, 0), ('', 18, 1)], None)
}

class TestName(object):
    
    def test_name(self):
        
        nameproc = lipyd.name.LipidNameProcessor()
        
        for name, result in swl_names.items():
            
            assert (
                nameproc.process(name, database = 'swisslipids') ==
                result
            )
        
        for name, result in lmp_names.items():
            
            assert (
                nameproc.process(name, database = 'lipidmaps') ==
                result
            )
