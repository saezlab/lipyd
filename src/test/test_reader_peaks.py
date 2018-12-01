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

import lipyd.reader.peaks
import lipyd.settings

class TestReaderPeaks(object):
    
    def test_reader_peaks(self):
        
        path = lipyd.settings.get('peaks_example')
        
        reader = lipyd.reader.peaks.PeaksReader(path)
        samples = reader.samples
        
        assert reader.mzs.shape == (443, 7)
        assert ('A', 10) in {s['label']['sample_id'] for s in samples}
        assert len({s['Normalized Area'] - s['m/z'] for s in samples}) == 1
        assert all(
            (
                samples[i    ]['label']['sample_id'] <=
                samples[i + 1]['label']['sample_id']
            )
            for i in range(len(samples) - 1)
        )
