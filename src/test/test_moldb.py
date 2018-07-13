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
    
    def test_lipidmaps(self):
        
        lm = lipyd.moldb.LipidMaps()
        
        gibberellin = list(
            lm.get_record('LMPR0104170034', typ = 'mainkey')
        )[0]
        
        assert gibberellin['annot']['COMMON_NAME'] == 'gibberellin A17'
