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
    
    nameproc = lipyd.name.LipidNameProcessor()
    
    def test_name_swisslipids(self):
        
        for name, result in swl_names.items():
            
            assert (
                self.nameproc.process(name, database = 'swisslipids') ==
                result
            )
    
    def test_name_lipidmaps(self):
        
        for name, result in lmp_names.items():
            
            assert (
                self.nameproc.process(name, database = 'lipidmaps') ==
                result
            )
    
    def test_name_iso(self):
        
        self.nameproc.iso = True
        self.nameproc.database = 'lipidmaps'
        
        assert (
            self.nameproc.process(
                'PE(16:0/18:1(9Z))-15-isoLG hydroxylactam'
            ) ==
            (
                'PE',
                ['', 34, 1],
                [('', 16, 0), ('', 18, 1)],
                [(('', 16, 0), '', False), (('', 18, 1), '9Z', False)]
            )
        )
        
        assert (
            self.nameproc.process('FAHFA(16:0/10-O-18:0)') ==
            (
                'FA',
                ['10-O', 34, 0],
                [('', 16, 0), ('10-O', 18, 0)],
                [(('', 16, 0), '', False), (('10-O', 18, 0), '', False)]
            )
        )
