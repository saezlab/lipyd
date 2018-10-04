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

import lipyd.reader.peaks
import lipyd.settings

class TestReaderPeaks(object):
    
    def test_reader_peaks(self):
        
        path = lipyd.settings.get('peaks_example')
        
        reader = lipyd.reader.peaks.PeaksReader(path)
        
        assert reader.mzs.shape == (443, 7)
