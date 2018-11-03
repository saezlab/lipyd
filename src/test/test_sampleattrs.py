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

import re

import lipyd.sampleattrs as sampleattrs
from lipyd.common import basestring

class TestSampleAttrs(object):
    
    def test_plate_sample_id_processor(self):
        
        sip = sampleattrs.plate_sample_id_processor()
        
        a11_from_string = sip('A11')
        a11_from_tuple  = sip(('A', 11))
        
        assert a11_from_string == ('A', 11)
        assert a11_from_tuple  == ('A', 11)
    
    def test_generic_sample_id_processor(self):
        
        def _method(something):
            
            if isinstance(something, basestring):
                
                retime = re.compile(r'([0-9]+)\s?(s|min|h)')
                m = retime.search(something)
                
                if m:
                    
                    m = m.groups()
                    
                    return (int(m[0]), m[1])
            
            return something, None

        sip = sampleattrs.sample_id_processor(_method, 'time', 'unit')

        min10 = sip('10 min')
        h6    = sip('6h')
        other = sip('foobar')
        
        assert min10 == (10, 'min')
        assert h6    == (6,  'h')
        assert other == ('foobar', None)
    
    def test_passthrough_sample_id_processor(self):
        
        sip = sampleattrs.sample_id_processor()
        
        foobar = sip('foobar')
        
        assert hasattr(foobar, 'sample_id')
        assert foobar.sample_id == 'foobar'
