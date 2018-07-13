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

import lipyd.moldb


class TestMoldb(object):
    
    lm = lipyd.moldb.LipidMaps()
    
    def test_lipidmaps_record(self):
        
        gibberellin = list(
            self.lm.get_record('LMPR0104170034', typ = 'mainkey')
        )[0]
        
        assert gibberellin['name']['COMMON_NAME'] == 'gibberellin A17'
    
    def test_lipidmaps_obmol(self):
        
        tag = list(self.lm.get_obmol('TAG(15:0_20:4_20:5)', 'synonym'))[0]
        
        assert abs(tag.exactmass - 886.7050407280012) < 0.000001
