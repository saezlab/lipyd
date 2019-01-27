#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `lipyd` python module
#
#  Copyright (c) 2015-2019 - EMBL
#
#  File author(s):
#  Dénes Türei (turei.denes@gmail.com)
#  Igor Bulanov
#
#  Distributed under the GNU GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://www.ebi.ac.uk/~denes
#

import pytest

import lipyd.settings as settings
import lipyd.sec as sec


class TestSec(object):
    """ """
    
    def test_sec_unicorn_asc(self):
        """ """
        
        path = settings.get('sec_unicorn_example')
        reader = sec.SECReader(path)
        
        highest = sorted(
            reader.profile(),
            key = lambda fr: fr.mean,
            reverse = True,
        )[0]
        
        assert highest.row == 'A' and highest.col == 10
    
    def test_sec_xls(self):
        """ """
        
        path = settings.get('sec_xls_example')
        reader = sec.SECReader(path)
        
        highest015 = sorted(
            reader.profile(start_volume = 0.615),
            key = lambda fr: fr.mean,
            reverse = True,
        )[0]
        
        highest045 = sorted(
            reader.profile(start_volume = 0.645),
            key = lambda fr: fr.mean,
            reverse = True,
        )[0]
        
        assert highest015.row == 'A' and highest015.col == 12
        assert highest045.row == 'A' and highest045.col == 11
