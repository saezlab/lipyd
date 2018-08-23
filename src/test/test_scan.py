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
                    score = 10,
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
    # PE positive
    (
        'pos_examples.mgf',
        'pos',
        2756,
        {'PE(34:1)': (
            MS2Identity(
                score = 15,
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
    # PC negative
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
    # Lyso-PE negative
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
    # Lyso-PE positive
    (
        'pos_examples.mgf',
        'pos',
        619,
        {'Lyso-PE(16:0)': (
            MS2Identity(
                score = 15,
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
    ),
    # Retinol positive
    (
        'pos_examples.mgf',
        'pos',
        1380,
        {'VA': (
            MS2Identity(
                score = 8,
                hg = Headgroup(main='VA', sub=()),
                chainsum = None,
                chains = None,
                details = None,
            ),
        )}
    ),
    # Retinol negative
    (
        'neg_examples.mgf',
        'neg',
        673,
        {'VA': (
            MS2Identity(
                score = 8,
                hg = Headgroup(main='VA', sub=()),
                chainsum = None,
                chains = None,
                details = None,
            ),
        )}
    ),
    # dCer positive
    (
        'pos_examples.mgf',
        'pos',
        3172,
        {'Cer(d36:1)': (
            MS2Identity(
                score = 19,
                hg = Headgroup(main='Cer', sub=()),
                chainsum = ChainSummary(
                    c = 36,
                    u = 1,
                    typ = ('Sph', 'FA'),
                    attr = (
                        ChainAttr(sph='d', ether=False, oh=()),
                        ChainAttr(sph='', ether=False, oh=())
                    )
                ),
                chains = (
                    Chain(
                        c = 18,
                        u = 1,
                        typ = 'Sph',
                        attr = ChainAttr(sph='d', ether=False, oh=()),
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
    # DHCer positive
    (
        'pos_examples.mgf',
        'pos',
        2367,
        {'Cer(DH36:0)': (
            MS2Identity(
                score = 22,
                hg = Headgroup(main='Cer', sub=()),
                chainsum = ChainSummary(
                    c = 36,
                    u = 0,
                    typ = ('Sph', 'FA'),
                    attr = (
                        ChainAttr(sph='DH', ether=False, oh=()),
                        ChainAttr(sph='', ether=False, oh=())
                    )
                ),
                chains = (
                    Chain(
                        c = 18,
                        u = 0,
                        typ = 'Sph',
                        attr = ChainAttr(sph='DH', ether=False, oh=()),
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
    # tCer positive
    (
        'pos_examples.mgf',
        'pos',
        2070,
        {'Cer(t36:0)': (
            MS2Identity(
                score = 25,
                hg = Headgroup(main='Cer', sub=()),
                chainsum = ChainSummary(
                    c = 36,
                    u = 0,
                    typ = ('Sph', 'FA'),
                    attr = (
                        ChainAttr(sph='t', ether=False, oh=()),
                        ChainAttr(sph='', ether=False, oh=())
                    )
                ),
                chains = (
                    Chain(
                        c = 18,
                        u = 0,
                        typ = 'Sph',
                        attr = ChainAttr(sph='t', ether=False, oh=()),
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
    # 2OH-acyl-dCer positive
    (
        'pos_examples.mgf',
        'pos',
        2421,
        {'Cer(d36:1-2OH)': (
            MS2Identity(
                score = 17,
                hg = Headgroup(main='Cer', sub=()),
                chainsum = ChainSummary(
                    c = 36,
                    u = 1,
                    typ = ('Sph', 'FA'),
                    attr = (
                        ChainAttr(sph='d', ether=False, oh=()),
                        ChainAttr(sph='', ether=False, oh=('2OH',))
                    )
                ),
                chains = (
                    Chain(
                        c = 18,
                        u = 1,
                        typ = 'Sph',
                        attr = ChainAttr(sph='d', ether=False, oh=()),
                        iso = ()
                    ),
                    Chain(
                        c = 18,
                        u = 0,
                        typ = 'FA',
                        attr = ChainAttr(sph='', ether=False, oh=('2OH',)),
                        iso = ()
                    )
                ),
                details = None
            ),
        )}
    ),
    # dCer-1-P positive
    # we have tCer-1-P neither in standards nor in screens hence
    # we can not have a test for that
    (
        'pos_examples.mgf',
        'pos',
        3208,
        {'Cer1P(d42:2)': (
            MS2Identity(
                score = 22,
                hg = Headgroup(main='Cer', sub=('1P',)),
                chainsum = ChainSummary(
                    c = 42,
                    u = 2,
                    typ = ('Sph', 'FA'),
                    attr = (
                        ChainAttr(sph='d', ether=False, oh=()),
                        ChainAttr(sph='', ether=False, oh=())
                    )
                ),
                chains = (
                    Chain(
                        c = 18,
                        u = 1,
                        typ = 'Sph',
                        attr = ChainAttr(sph='d', ether=False, oh=()),
                        iso = ()
                    ),
                    Chain(
                        c = 24,
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
    # Hex-dCer positive
    (
        'pos_examples.mgf',
        'pos',
        3806,
        {'Hex-Cer(d42:2)': (
            MS2Identity(
                score = 35,
                hg = Headgroup(main='Cer', sub=('Hex',)),
                chainsum = ChainSummary(
                    c = 42,
                    u = 2,
                    typ = ('Sph', 'FA'),
                    attr = (
                        ChainAttr(sph='d', ether=False, oh=()),
                        ChainAttr(sph='', ether=False, oh=())
                    )
                ),
                chains = (
                    Chain(
                        c = 18,
                        u = 1,
                        typ = 'Sph',
                        attr = ChainAttr(sph='d', ether=False, oh=()),
                        iso = ()
                    ),
                    Chain(
                        c = 24,
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
    # Hex-tCer positive
    (
        'pos_examples.mgf',
        'pos',
        3787,
        {'Hex-Cer(t42:2)': (
            MS2Identity(
                score = 36,
                hg = Headgroup(main='Cer', sub=('Hex',)),
                chainsum = ChainSummary(
                    c = 42,
                    u = 2,
                    typ = ('Sph', 'FA'),
                    attr = (
                        ChainAttr(sph='t', ether=False, oh=()),
                        ChainAttr(sph='', ether=False, oh=())
                    )
                ),
                chains = (
                    Chain(
                        c = 18,
                        u = 1,
                        typ = 'Sph',
                        attr = ChainAttr(sph='t', ether=False, oh=()),
                        iso = ()
                    ),
                    Chain(
                        c = 24,
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
    # Hex2-tCer positive
    (
        'pos_examples.mgf',
        'pos',
        3626,
        {'Hex2-Cer(t42:2)': (
            MS2Identity(
                score = 42,
                hg = Headgroup(main='Cer', sub=('Hex2',)),
                chainsum = ChainSummary(
                    c = 42,
                    u = 2,
                    typ = ('Sph', 'FA'),
                    attr = (
                        ChainAttr(sph='t', ether=False, oh=()),
                        ChainAttr(sph='', ether=False, oh=())
                    )
                ),
                chains = (
                    Chain(
                        c = 18,
                        u = 1,
                        typ = 'Sph',
                        attr = ChainAttr(sph='t', ether=False, oh=()),
                        iso = ()
                    ),
                    Chain(
                        c = 24,
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
    # Hex-tCer positive
    (
        'pos_examples.mgf',
        'pos',
        3005,
        {'Hex-Cer(t40:1)': (
            MS2Identity(
                score = 26,
                hg = Headgroup(main='Cer', sub=('Hex',)),
                chainsum = ChainSummary(
                    c = 40,
                    u = 1,
                    typ = ('Sph', 'FA'),
                    attr = (
                        ChainAttr(sph='t', ether=False, oh=()),
                        ChainAttr(sph='', ether=False, oh=())
                    )
                ),
                chains = (
                    Chain(
                        c = 18,
                        u = 0,
                        typ = 'Sph',
                        attr = ChainAttr(sph='t', ether=False, oh=()),
                        iso = ()
                    ),
                    Chain(
                        c = 22,
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
    # Hex-tCer positive
    (
        'pos_examples.mgf',
        'pos',
        2981,
        {'Hex-Cer(t42:2)': (
            MS2Identity(
                score = 26,
                hg = Headgroup(main='Cer', sub=('Hex',)),
                chainsum = ChainSummary(
                    c = 42,
                    u = 2,
                    typ = ('Sph', 'FA'),
                    attr = (
                        ChainAttr(sph='t', ether=False, oh=()),
                        ChainAttr(sph='', ether=False, oh=())
                    )
                ),
                chains = (
                    Chain(
                        c = 18,
                        u = 0,
                        typ = 'Sph',
                        attr = ChainAttr(sph='t', ether=False, oh=()),
                        iso = ()
                    ),
                    Chain(
                        c = 24,
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
    # Hex-dCer positive
    (
        'pos_examples.mgf',
        'pos',
        3092,
        {'Hex-Cer(d42:2)': (
            MS2Identity(
                score = 25,
                hg = Headgroup(main='Cer', sub=('Hex',)),
                chainsum = ChainSummary(
                    c = 42,
                    u = 2,
                    typ = ('Sph', 'FA'),
                    attr = (
                        ChainAttr(sph='d', ether=False, oh=()),
                        ChainAttr(sph='', ether=False, oh=())
                    )
                ),
                chains = (
                    Chain(
                        c = 18,
                        u = 1,
                        typ = 'Sph',
                        attr = ChainAttr(sph='d', ether=False, oh=()),
                        iso = ()
                    ),
                    Chain(
                        c = 24,
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
        2290,
        {'SM(d34:1)': (
            MS2Identity(
                score = 35,
                hg = Headgroup(main='SM', sub=()),
                chainsum = ChainSummary(
                    c = 34,
                    u = 1,
                    typ = ('Sph', 'FA'),
                    attr = (
                        ChainAttr(sph='d', ether=False, oh=()),
                        ChainAttr(sph='', ether=False, oh=())
                    )
                ),
                chains = (
                    Chain(
                        c = 18,
                        u = 1,
                        typ = 'Sph',
                        attr = ChainAttr(sph='d', ether=False, oh=()),
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
    # SHex-dCer
    (
        'pos_examples.mgf',
        'pos',
        2550,
        {'SHex-Cer(d42:2)': (
            MS2Identity(
                score = 26,
                hg = Headgroup(main='Cer', sub=('SHex',)),
                chainsum = ChainSummary(
                    c = 42,
                    u = 2,
                    typ = ('Sph', 'FA'),
                    attr = (
                        ChainAttr(sph='d', ether=False, oh=()),
                        ChainAttr(sph='', ether=False, oh=())
                    )
                ),
                chains = (
                    Chain(
                        c = 18,
                        u = 1,
                        typ = 'Sph',
                        attr = ChainAttr(sph='d', ether=False, oh=()),
                        iso = ()
                    ),
                    Chain(
                        c = 24,
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
    # dSph positive
    (
        'pos_examples.mgf',
        'pos',
        566,
        {'Sph(d18:1)': (
            MS2Identity(
                score = 27,
                hg = Headgroup(main='Sph', sub=()),
                chainsum = ChainSummary(
                    c = 18,
                    u = 1,
                    typ = ('Sph',),
                    attr = (ChainAttr(sph='d', ether=False, oh=()),)
                ),
                chains = (
                    Chain(
                        c = 18,
                        u = 1,
                        typ = 'Sph',
                        attr = ChainAttr(sph='d', ether=False, oh=()),
                        iso = ()
                    ),
                ),
                details = None
            ),
        )}
    ),
    # DHSph positive
    (
        'pos_examples.mgf',
        'pos',
        453,
        {'Sph(DH18:0)': (
            MS2Identity(
                score = 27,
                hg = Headgroup(main='Sph', sub=()),
                chainsum = ChainSummary(
                    c = 18,
                    u = 0,
                    typ = ('Sph',),
                    attr = (ChainAttr(sph='DH', ether=False, oh=()),)
                ),
                chains = (
                    Chain(
                        c = 18,
                        u = 0,
                        typ = 'Sph',
                        attr = ChainAttr(sph='DH', ether=False, oh=()),
                        iso = ()
                    ),
                ),
                details = None
            ),
        )}
    ),
    # tSph positive
    (
        'pos_examples.mgf',
        'pos',
        354,
        {'Sph(t18:0)': (
            MS2Identity(
                score = 33,
                hg = Headgroup(main='Sph', sub=()),
                chainsum = ChainSummary(
                    c = 18,
                    u = 0,
                    typ = ('Sph',),
                    attr = (ChainAttr(sph='t', ether=False, oh=()),)
                ),
                chains = (
                    Chain(
                        c = 18,
                        u = 0,
                        typ = 'Sph',
                        attr = ChainAttr(sph='t', ether=False, oh=()),
                        iso = ()
                    ),
                ),
                details = None
            ),
        )}
    ),
    # kSph positive
    (
        'pos_examples.mgf',
        'pos',
        48,
        {'Sph(k8:0)': (
            MS2Identity(
                score = 50,
                hg = Headgroup(main='Sph', sub=()),
                chainsum = ChainSummary(
                    c = 8,
                    u = 0,
                    typ = ('Sph',),
                    attr = (ChainAttr(sph='k', ether=False, oh=()),)
                ),
                chains = (
                    Chain(
                        c = 8,
                        u = 0,
                        typ = 'Sph',
                        attr = ChainAttr(sph='k', ether=False, oh=()),
                        iso = ()
                    ),
                ),
                details = None
            ),
        )}
    ),
    # DHSph-1P positive
    (
        'pos_examples.mgf',
        'pos',
        465,
        {'': (

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
