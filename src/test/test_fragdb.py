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

import os
import numpy as np

import lipyd.mgf as mgf
import lipyd.fragdb as fragdb
import lipyd.ms2 as ms2
import lipyd.settings as settings


class TestFragdb(object):
    
    mgffile = settings.get('mgf_pos_examples')
    mgfreader = mgf.MgfReader(mgffile)
    
    def test_mgf_reader(self):
        
        precursor = 590.45536 # this is a Cer-1P
        idx, rtdiff = self.mgfreader.lookup(precursor)
        
        tolerance = (
            precursor / 1e06 * settings.get('precursor_match_tolerance')
        )
        
        assert np.all(
            np.abs(
                self.mgfreader.mgfindex[idx,0] - precursor
            ) <= tolerance
        )
    
    def test_annotate(self):
        
        precursor = 590.45536 # this is a Cer-1P(32:1)
        scan = self.mgfreader.scan_by_id(1941)
        
        annot = fragdb.FragmentAnnotator(
            mzs = scan[:,0],
            ionmode = 'pos',
            precursor = precursor
        )
        
        fragnames = set(aa.name for a in annot for aa in a)
        
        assert '[FA(14:0)+NH+C2H2-OH]+' in fragnames
        assert '[Sph(18:1)-2xH2O+H]+' in fragnames
        assert len(list(annot)) == len(annot.mzs)
