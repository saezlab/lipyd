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

import os

import lipyd.mgf as mgf
import lipyd.fragdb as fragdb
import lipyd.ms2 as ms2
import lipyd.settings as settings
import lipyd.common as common
import lipyd.moldb as moldb
import lipyd.lipproc as lipproc
from lipyd.lipproc import Headgroup, Chain, ChainSummary, ChainAttr
from lipyd.ms2 import MS2Identity

# do not include details now it would be too much
settings.setup(ms2_scan_chain_details = False)


specimens = [
    # FA negative mode
    (
        'neg_examples.mgf',
        'neg',
        691,
        {
            'FA(18:1)': (
                ms2.MS2Identity(
                    score = 5,
                    hg = lipproc.Headgroup(main = 'FA'),
                    chainsum = lipproc.ChainSummary(
                        c = 18,
                        u = 1,
                        typ = ('FA',),
                        attr = (
                            lipproc.ChainAttr(
                                sph = '', ether = False, oh = ()
                            ),
                        )
                    ),
                    chains = (
                        lipproc.Chain(
                            c = 18,
                            u = 1,
                            typ = 'FA',
                            attr = lipproc.ChainAttr(
                                sph = '', ether = False, oh = ()
                            ),
                            iso = ()
                        ),
                    )
                ),
            )
        }
    ),
    (
        'neg_examples.mgf',
        'neg',
        134,
        {}
    ),
    # DAG positive, in vivo SEC14L2
    (
        'pos_examples.mgf',
        'pos',
        3344,
        {'DAG(32:1)': (
            MS2Identity(
                score = 6,
                hg = Headgroup(main='DAG', sub=()),
                chainsum = ChainSummary(
                    c = 32,
                    u = 1,
                    typ = ('FA', 'FA'),
                    attr = (
                        ChainAttr(sph='', ether=False, oh=()),
                        ChainAttr(sph='', ether=False, oh=())
                    )
                ),
                chains = (
                    Chain(
                        c = 16,
                        u = 0,
                        typ = 'FA',
                        attr = ChainAttr(sph='', ether=False, oh=()),
                        iso = ()
                    ),
                    Chain(
                        c = 16,
                        u = 1,
                        typ = 'FA',
                        attr = ChainAttr(sph='', ether=False, oh=()),
                        iso = ()
                    )
                ),
                details = None
            ),
        )}
    ),
    (
        'pos_examples.mgf',
        'pos',
        4923,
        {'TAG(48:3)': (
            MS2Identity(
                score = 5,
                hg = Headgroup(main='TAG', sub=()),
                chainsum = ChainSummary(
                    c = 48,
                    u = 3,
                    typ = ('FA', 'FA', 'FA'),
                    attr = (
                        ChainAttr(sph='', ether=False, oh=()),
                        ChainAttr(sph='', ether=False, oh=()),
                        ChainAttr(sph='', ether=False, oh=())
                    ),
                    iso = None
                ),
                chains = (
                    Chain(
                        c = 16,
                        u = 1,
                        typ = 'FA',
                        attr = ChainAttr(sph='', ether=False, oh=()),
                        iso = ()
                    ),
                    Chain(
                        c = 16,
                        u = 1,
                        typ = 'FA',
                        attr = ChainAttr(sph='', ether=False, oh=()),
                        iso = ()
                    ),
                    Chain(
                        c = 16,
                        u = 1,
                        typ = 'FA',
                        attr = ChainAttr(sph='', ether=False, oh=()),
                        iso = ()
                    )
                ),
                details = None
            ),
            MS2Identity(
                score = 5,
                hg = Headgroup(main='TAG', sub=()),
                chainsum = ChainSummary(
                    c = 48,
                    u = 3,
                    typ = ('FA', 'FA', 'FA'),
                    attr = (
                        ChainAttr(sph='', ether=False, oh=()),
                        ChainAttr(sph='', ether=False, oh=()),
                        ChainAttr(sph='', ether=False, oh=())
                    ),
                    iso = None
                ),
                chains = (
                    Chain(
                        c = 18,
                        u = 1,
                        typ = 'FA',
                        attr = ChainAttr(sph='', ether=False, oh=()),
                        iso = ()
                    ),
                    Chain(
                        c = 14,
                        u = 1,
                        typ = 'FA',
                        attr = ChainAttr(sph='', ether=False, oh=()),
                        iso = ()
                    ),
                    Chain(
                        c = 16,
                        u = 1,
                        typ = 'FA',
                        attr = ChainAttr(sph='', ether=False, oh=()),
                        iso = ()
                    )
                ),
                details = None
            )
        )}
    )
]


class TestScan(object):
    
    @pytest.mark.parametrize('mgfname, ionmode, scan_id, identity', specimens)
    def test_identification(
            self,
            mgfname,
            ionmode,
            scan_id,
            identity
        ):
        
        mgfpath = os.path.join(common.ROOT, 'data', 'ms2_examples', mgfname)
        scan = ms2.Scan.from_mgf(
            mgfpath, scan_id, ionmode
        )
        
        iscan = scan.identify()
        
        for name, ids in identity.items():
            
            assert name in iscan
            
            for i in ids:
                
                assert i in iscan[name]
