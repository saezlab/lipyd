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
    # TAG positive, in vivo STARD11 + 818.7187
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
                    )
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
                    )
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
    ),
    (
        'neg_examples.mgf',
        'neg',
        2961,
        {'PE(34:1)': (
            MS2Identity(
                score = 6.5,
                hg = Headgroup(main='PE', sub=()),
                chainsum = ChainSummary(
                    c = 34,
                    u = 1,
                    typ = ('FA', 'FA'),
                    attr = (
                        ChainAttr(sph='', ether=False, oh=()),
                        ChainAttr(sph='', ether=False, oh=())
                    )
                ),
                chains = (
                    Chain(
                        c = 18,
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
        2756,
        {'PE(34:1)': (
            MS2Identity(
                score = 5,
                hg = Headgroup(main='PE', sub=()),
                chainsum = ChainSummary(
                    c = 34,
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
                        c = 18,
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
        'neg_examples.mgf',
        'neg',
        2329,
        {'PC(32:1)': (
            MS2Identity(
                score = 6.0,
                hg = Headgroup(main='PC', sub=()),
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
        ),
        # this to see if score is 0
        'PE(34:1)': (
            MS2Identity(
                score = 0,
                hg = Headgroup(main='PE', sub=()),
                chainsum = ChainSummary(
                    c = 34,
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
                        c = 18,
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
        2731,
        {'PC(36:2)': (
            MS2Identity(
                score = 9,
                hg = Headgroup(main='PC', sub=()),
                chainsum = ChainSummary(
                    c = 36,
                    u = 2,
                    typ = ('FA', 'FA'),
                    attr = (
                        ChainAttr(sph='', ether=False, oh=()),
                        ChainAttr(sph='', ether=False, oh=())
                    )
                ),
                chains = None,
                details = None
            ),
        )}
    ),
    (
        'neg_examples.mgf',
        'neg',
        2516,
        {'PI(34:1)': (
            MS2Identity(
                score = 11.0,
                hg = Headgroup(main='PI', sub=()),
                chainsum = ChainSummary(
                    c = 34,
                    u = 1,
                    typ = ('FA', 'FA'),
                    attr = (
                        ChainAttr(sph='', ether=False, oh=()),
                        ChainAttr(sph='', ether=False, oh=())
                    )
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
                        c = 16,
                        u = 0,
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
        2917,
        {'PI(38:3)': (
            MS2Identity(
                score = 5,
                hg = Headgroup(main='PI', sub=()),
                chainsum = ChainSummary(
                    c = 38,
                    u = 3,
                    typ = ('FA', 'FA'),
                    attr = (
                        ChainAttr(sph='', ether=False, oh=()),
                        ChainAttr(sph='', ether=False, oh=())
                    )
                ),
                chains = (
                    Chain(
                        c = 20,
                        u = 3,
                        typ = 'FA',
                        attr = ChainAttr(sph='', ether=False, oh=()),
                        iso = ()
                    ),
                    Chain(
                        c = 18,
                        u = 0,
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
        'neg_examples.mgf',
        'neg',
        2738,
        {'PS(36:1)': (
            MS2Identity(
                score = 7.0,
                hg = Headgroup(main='PS', sub=()),
                chainsum = ChainSummary(
                    c = 36,
                    u = 1,
                    typ = ('FA', 'FA'),
                    attr = (
                        ChainAttr(sph='', ether=False, oh=()),
                        ChainAttr(sph='', ether=False, oh=())
                    )
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
                        c = 18,
                        u = 0,
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
        2759,
        {'PS(36:1)': (
            MS2Identity(
                score = 5,
                hg = Headgroup(main='PS', sub=()),
                chainsum = ChainSummary(
                    c = 36,
                    u = 1,
                    typ = ('FA', 'FA'),
                    attr = (
                        ChainAttr(sph='', ether=False, oh=()),
                        ChainAttr(sph='', ether=False, oh=())
                    )
                ),
                chains = (
                    Chain(
                        c = 18,
                        u = 0,
                        typ = 'FA',
                        attr = ChainAttr(sph='', ether=False, oh=()),
                        iso = ()
                    ),
                    Chain(
                        c = 18,
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
        'neg_examples.mgf',
        'neg',
        2551,
        {'PG(38:3)': (
            MS2Identity(
                score = 5.0,
                hg = Headgroup(main='PG', sub=()),
                chainsum = ChainSummary(
                    c = 38,
                    u = 3,
                    typ = ('FA', 'FA'),
                    attr = (
                        ChainAttr(sph='', ether=False, oh=()),
                        ChainAttr(sph='', ether=False, oh=())
                    )
                ),
                chains = (
                    Chain(
                        c = 20,
                        u = 2,
                        typ = 'FA',
                        attr = ChainAttr(sph='', ether=False, oh=()),
                        iso = ()
                    ),
                    Chain(
                        c = 18,
                        u = 1,
                        typ = 'FA',
                        attr = ChainAttr(sph='', ether=False, oh=()),
                        iso = ()
                    )
                ),
                details = None
            ),
        ),
        'BMP(38:3)': (
            MS2Identity(
                score = 5.0,
                hg = Headgroup(main='BMP', sub=()),
                chainsum = ChainSummary(
                    c = 38,
                    u = 3,
                    typ = ('FA', 'FA'),
                    attr = (
                        ChainAttr(sph='', ether=False, oh=()),
                        ChainAttr(sph='', ether=False, oh=())
                    )
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
                        c = 20,
                        u = 2,
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
        2552,
        {'BMP(36:2)': (
            MS2Identity(
                score = 5,
                hg = Headgroup(main='BMP', sub=()),
                chainsum = ChainSummary(
                    c = 36,
                    u = 2,
                    typ = ('FA', 'FA'),
                    attr = (
                        ChainAttr(sph='', ether=False, oh=()),
                        ChainAttr(sph='', ether=False, oh=())
                    )
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
                        c = 18,
                        u = 1,
                        typ = 'FA',
                        attr = ChainAttr(sph='', ether=False, oh=()),
                        iso = ()
                    )
                ),
                details = None
            ),
        ),
        'PG(36:2)': (
            MS2Identity(
                score = 0,
                hg = Headgroup(main='PG', sub=()),
                chainsum = ChainSummary(
                    c = 36,
                    u = 2,
                    typ = ('FA', 'FA'),
                    attr = (
                        ChainAttr(sph='', ether=False, oh=()),
                        ChainAttr(sph='', ether=False, oh=())
                    )
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
                        c = 18,
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
    # Lyso-PE
    (
        'neg_examples.mgf',
        'neg',
        361,
        {'Lyso-PE(16:1)': (
            MS2Identity(
                score = 6.0,
                hg = Headgroup(main='PE', sub=('Lyso',)),
                chainsum = ChainSummary(
                    c = 16,
                    u = 1,
                    typ = ('FA',),
                    attr = (ChainAttr(sph='', ether=False, oh=()),)
                ),
                chains = (
                    Chain(
                        c = 16,
                        u = 1,
                        typ = 'FA',
                        attr = ChainAttr(sph='', ether=False, oh=()),
                        iso = ()
                    ),
                ),
                details = None
            ),
        )}
    ),
    (
        'pos_examples.mgf',
        'pos',
        619,
        {'Lyso-PE(16:0)': (
            MS2Identity(
                score = 10,
                hg = Headgroup(main='PE', sub=('Lyso',)),
                chainsum = ChainSummary(
                    c = 16,
                    u = 0,
                    typ = ('FA',),
                    attr = (ChainAttr(sph='', ether=False, oh=()),)
                ),
                chains = (
                    Chain(
                        c = 16,
                        u = 0,
                        typ = 'FA',
                        attr = ChainAttr(sph='', ether=False, oh=()),
                        iso = ()
                    ),
                ),
                details = None
            ),
        )}
    ),
    # Lyso-PC
    (
        'pos_examples.mgf',
        'pos',
        604,
        {'Lyso-PC(18:1)': (
            MS2Identity(
                score = 10,
                hg = Headgroup(main='PC', sub=('Lyso',)),
                chainsum = ChainSummary(
                    c = 18,
                    u = 1,
                    typ = ('FA',),
                    attr = (ChainAttr(sph='', ether=False, oh=()),)
                ),
                chains = (
                    Chain(
                        c = 18,
                        u = 1,
                        typ = 'FA',
                        attr = ChainAttr(sph='', ether=False, oh=()),
                        iso = ()
                    ),
                ),
                details = None
            ),
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
        scan = ms2.Scan.from_mgf(mgfpath, scan_id, ionmode)
        
        iscan = scan.identify()
        
        for name, ids in identity.items():
            
            assert name in iscan
            
            for i in ids:
                
                assert i in iscan[name]
