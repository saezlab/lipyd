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

import lipyd.name as name
import lipyd.lipproc as lipproc


class TestLipproc(object):
    
    def test_is_subset_of(self):
        
        empty_lab = lipproc.LipidLabel(None, None, None, None)
        lnp = name.LipidNameProcessor()
        
        pe = lnp.process('Phosphatidylethanolamine')
        pe_361 = lnp.process('Phosphatidylethanolamine(36:1)')
        pe_180181 = lnp.process('Phosphatidylethanolamine(18:0/18:1)')
        pe = lipproc.LipidRecord(
            lab = empty_lab,
            hg = pe[0],
            chainsum = None,
            chains = (),
        )
        pe_361 = lipproc.LipidRecord(
            lab = empty_lab,
            hg = pe_361[0],
            chainsum = pe_361[1],
            chains = (),
        )
        pe_180181 = lipproc.LipidRecord(
            lab = empty_lab,
            hg = pe_180181[0],
            chainsum = pe_180181[1],
            chains = pe_180181[2],
        )
        
        assert pe_361 in pe
        assert pe_180181 in pe_361
        assert pe_180181 in pe
        
        cer_d_361 = lnp.process('Phosphatidylethanolamine(d36:1)')
        cer_dh_361 = lnp.process('Phosphatidylethanolamine(DH36:1)')
        cer_d_361 = lipproc.LipidRecord(
            lab = empty_lab,
            hg = cer_d_361[0],
            chainsum = cer_d_361[1],
            chains = (),
        )
        cer_dh_361 = lipproc.LipidRecord(
            lab = empty_lab,
            hg = cer_dh_361[0],
            chainsum = cer_dh_361[1],
            chains = (),
        )
        
        assert cer_d_361 not in cer_dh_361
        assert cer_dh_361 not in cer_d_361
