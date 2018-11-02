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

import itertools
import numpy as np
import lipyd.moldb
import lipyd.lipproc as lipproc


class TestMoldb(object):
    
    lm = lipyd.moldb.LipidMaps()
    mda = lipyd.moldb.MoleculeDatabaseAggregator()
    
    def test_database_size(self):
        
        assert self.mda.data.shape[0] > 133000
    
    def test_chainlength_unsaturation(self):
        
        supposed_to_have = set(
            cu for cu in
            itertools.product(range(8, 73), range(0, 21))
            if cu[1] <= cu[0] - 6
        )
        
        actually_got = set(
            (r.chainsum.c, r.chainsum.u)
            for r in self.mda.data
            if (
                r.lab.db == 'lipyd.lipid' and
                r.hg.main == 'PC' and
                not r.hg.sub
            )
        )
        
        assert supposed_to_have == actually_got
    
    def test_lipidmaps_record(self):
        
        gibberellin = list(
            self.lm.get_record('LMPR0104170034', typ = 'mainkey')
        )[0]
        
        assert gibberellin['name']['COMMON_NAME'] == 'gibberellin A17'
    
    def test_lipidmaps_obmol(self):
        
        tag = list(self.lm.get_obmol('TAG(15:0_20:4_20:5)', 'synonym'))[0]
        
        assert abs(tag.exactmass - 886.7050407280012) < 0.000001
    
    def test_aggregator_build(self):
        
        assert len([i for i in self.mda.masses if i == 0 or np.isnan(i)]) == 0
        assert self.mda.masses.shape[0] == self.mda.data.shape[0]
    
    def test_aggregator_lookup(self):
        
        result = self.mda.lookup(808.634583)
        
        assert ('lipyd.lipid', 'PA', 45, 4) in [
            (i.lab.db, i.hg.main, i.chainsum.c, i.chainsum.u)
            for i in result[1]
            if i.hg is not None and i.chainsum is not None
        ]
        
        assert 'SLM:000056288' in [i.lab.db_id for i in result[1]]
    
    def test_aggregator_adduct_lookup(self):
        
        result = self.mda.adduct_lookup(
            728.605042778354, ionmode = 'pos'
        )['[M+H]+'][1]
        
        assert 'SLM:000391523' in [i.lab.db_id for i in result]
        
        swl_cer1p = [i for i in result if i.lab.db_id == 'SLM:000391523'][0]
        
        assert 'Ceramide phosphate (d42:2)' in swl_cer1p.lab.names
        assert '' not in swl_cer1p.lab.names
        assert '1P' in swl_cer1p.hg.sub
        
        lyp_cer1p = lipproc.LipidRecord(
            lab = lipproc.LipidLabel(
                db_id = None,
                db = 'lipyd.lipid',
                names = ('Cer-1P(DH42:2)',)
            ),
            hg = lipproc.Headgroup(
                main = 'Cer',
                sub = ('1P',)
            ),
            chainsum = lipproc.ChainSummary(
                c = 42, u = 2, typ = ('Sph', 'FA'),
                attr = (
                    lipproc.ChainAttr(sph = 'DH', ether = False, oh = ()),
                    lipproc.ChainAttr(sph = '', ether = False, oh = ())
                ),
            ),
            chains = ()
        )
        
        assert lyp_cer1p in list(result)
