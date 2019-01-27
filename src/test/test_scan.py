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


# certain cases for example when our standard suffers from a
# mass error, this way we apply special settings, e.g. increasing
# the tolerance at fragment mass lookups
special_settings = {
    # finally this was not necessary
    #(465, 'pos'): {
        #'scan_args': {'tolerance = 250'},
    #},
}


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
            )
        )}
    ),
    # PE negative, standards
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
            ),
        )}
    ),
    # PE positive, standards
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
            ),
        ),
        # this to see if score is 0
        'PE(34:1)': (
            False,
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
            ),
        ),
        'PG(36:2)': (
            False,
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
                    typ = ('Sph', 'FAOH'),
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
                        typ = 'FAOH',
                        attr = ChainAttr(sph='', ether=False, oh=('2OH',)),
                        iso = ()
                    )
                ),
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
        {'Cer-1P(d42:2)': (
            MS2Identity(
                score = 32,
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
                score = 39,
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
                score = 45,
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
            ),
        )}
    ),
    # Hex-tCer positive
    (
        'pos_examples.mgf',
        'pos',
        3005,
        {'Hex-Cer(d40:1-2OH)': (
            MS2Identity(
                score = 26,
                hg = Headgroup(main='Cer', sub=('Hex',)),
                chainsum = ChainSummary(
                    c = 40,
                    u = 1,
                    typ = ('Sph', 'FAOH'),
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
                        c = 22,
                        u = 0,
                        typ = 'FAOH',
                        attr = ChainAttr(sph='', ether=False, oh=('2OH',)),
                        iso = ()
                    )
                ),
            ),
        )}
    ),
    # Hex-tCer positive
    (
        'pos_examples.mgf',
        'pos',
        2981,
        {'Hex-Cer(d42:2-2OH)': (
            MS2Identity(
                score = 26,
                hg = Headgroup(main='Cer', sub=('Hex',)),
                chainsum = ChainSummary(
                    c = 42,
                    u = 2,
                    typ = ('Sph', 'FAOH'),
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
                        c = 24,
                        u = 1,
                        typ = 'FAOH',
                        attr = ChainAttr(sph='', ether=False, oh=('2OH',)),
                        iso = ()
                    )
                ),
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
            ),
        )}
    ),
    (
        'pos_examples.mgf',
        'pos',
        2290,
        {'SM(d34:1)': (
            MS2Identity(
                score = 41,
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
            ),
        )}
    ),
    # DHSph-1P positive
    # +180 ppm off
    (
        'pos_examples.mgf',
        'pos',
        568,
        {'Sph-1P(DH18:0)': (
            MS2Identity(
                score = 26,
                hg = Headgroup(main='Sph', sub=('1P',)),
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
            ),
        )}
    ),
    # tSph-1P positive
    (
        'pos_examples.mgf',
        'pos',
        338,
        {'Sph-1P(t18:0)': (
            MS2Identity(
                score = 38,
                hg = Headgroup(main='Sph', sub=('1P',)),
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
            ),
        )}
    ),
    # dSph-M2 positive
    (
        'pos_examples.mgf',
        'pos',
        473,
        {'Sph-M2(d18:1)': (
            MS2Identity(
                score = 65,
                hg = Headgroup(main='Sph', sub=('M2',)),
                chainsum = ChainSummary(
                    c = 18,
                    u = 1,
                    typ = ('Sph',),
                    attr = (ChainAttr(sph='d', ether=False, oh=()),),
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
            ),
        ),
        # to see if this gets lower score
        'Sph-M1(d19:1)': (
            False,
            MS2Identity(
                score = 32,
                hg = Headgroup(main='Sph', sub=('M1',)),
                chainsum = ChainSummary(
                    c = 19,
                    u = 1,
                    typ = ('Sph',),
                    attr = (ChainAttr(sph='d', ether=False, oh=()),),
                ),
                chains = (
                    Chain(
                        c = 19,
                        u = 1,
                        typ = 'Sph',
                        attr = ChainAttr(sph='d', ether=False, oh=()),
                        iso = ()
                    ),
                ),
            ),
        )},
    ),
    # dCer example from 2015, wrong collision energy
    (
        'neg_examples.mgf',
        'neg',
        3123,
        {'Cer(d36:1)': (
            MS2Identity(
                score = 51,
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
            ),
        ),
        'Cer(DH36:1)': (
            False,
            MS2Identity(
                score = 11,
                hg = Headgroup(main='Cer', sub=()),
                chainsum = ChainSummary(
                    c = 36,
                    u = 1,
                    typ = ('Sph', 'FA'),
                    attr = (
                        ChainAttr(sph='DH', ether=False, oh=()),
                        ChainAttr(sph='', ether=False, oh=())
                    )
                ),
                chains = (
                    Chain(
                        c = 18,
                        u = 1,
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
            ),
        )}
    ),
    # [DHCer-H]-
    (
        'neg_examples.mgf',
        'neg',
        2160,
        {'Cer(DH36:0)': (
            MS2Identity(
                score = 37,
                hg = Headgroup(main='Cer', sub=()),
                chainsum = ChainSummary(
                    c = 36,
                    u = 0,
                    typ = ('Sph', 'FA'),
                    attr = (
                        ChainAttr(sph='DH', ether=False, oh=()),
                        ChainAttr(sph='', ether=False, oh=())
                    ),
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
            ),
        )}
    ),
    # [DHCer+HCOO]-
    (
        'neg_examples.mgf',
        'neg',
        2153,
        {'Cer(DH36:0)': (
            MS2Identity(
                score = 34,
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
            ),
        )}
    ),
    # STARD11 negative dCer / DHCer formiate adducts
    (
        'neg_examples.mgf',
        'neg',
        2369,
        {'Cer(d32:1)': (
            MS2Identity(
                score = 54,
                hg = Headgroup(main='Cer', sub=()),
                chainsum = ChainSummary(
                    c = 32,
                    u = 1,
                    typ = ('Sph', 'FA'),
                    attr = (
                        ChainAttr(sph='d', ether=False, oh=()),
                        ChainAttr(sph='', ether=False, oh=())
                    ),
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
                        c = 14,
                        u = 0,
                        typ = 'FA',
                        attr = ChainAttr(sph='', ether=False, oh=()),
                        iso = ()
                    )
                ),
            ),
        )}
    ),
    (
        'neg_examples.mgf',
        'neg',
        2723,
        {'Cer(d34:1)': (
            MS2Identity(
                score = 51,
                hg = Headgroup(main='Cer', sub=()),
                chainsum = ChainSummary(
                    c = 34,
                    u = 1,
                    typ = ('Sph', 'FA'),
                    attr = (
                        ChainAttr(sph='d', ether=False, oh=()),
                        ChainAttr(sph='', ether=False, oh=())
                    ),
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
            ),
        )}
    ),
    # SEC14L1 negative dCer / DHCer formiate adducts
    (
        'neg_examples.mgf',
        'neg',
        2776,
        {'Cer(d34:1)': (
            MS2Identity(
                score = 54,
                hg = Headgroup(main='Cer', sub=()),
                chainsum = ChainSummary(
                    c = 34,
                    u = 1,
                    typ = ('Sph', 'FA'),
                    attr = (
                        ChainAttr(sph='d', ether=False, oh=()),
                        ChainAttr(sph='', ether=False, oh=())
                    ),
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
            ),
        )}
    ),
    (
        'neg_examples.mgf',
        'neg',
        2449,
        {'Cer(d32:1)': (
            MS2Identity(
                score = 48,
                hg = Headgroup(main='Cer', sub=()),
                chainsum = ChainSummary(
                    c = 32,
                    u = 1,
                    typ = ('Sph', 'FA'),
                    attr = (
                        ChainAttr(sph='d', ether=False, oh=()),
                        ChainAttr(sph='', ether=False, oh=())
                    ),
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
                        c = 14,
                        u = 0,
                        typ = 'FA',
                        attr = ChainAttr(sph='', ether=False, oh=()),
                        iso = ()
                    )
                ),
            ),
        )}
    ),
    (
        'neg_examples.mgf',
        'neg',
        2486,
        {'Cer(d34:2)': (
            MS2Identity(
                score = 31,
                hg = Headgroup(main='Cer', sub=()),
                chainsum = ChainSummary(
                    c = 34,
                    u = 2,
                    typ = ('Sph', 'FA'),
                    attr = (
                        ChainAttr(sph='d', ether=False, oh=()),
                        ChainAttr(sph='', ether=False, oh=())
                    )
                ),
                chains = (
                    Chain(
                        c = 16,
                        u = 1,
                        typ = 'Sph',
                        attr = ChainAttr(sph='d', ether=False, oh=()),
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
            ),
        )}
    ),
    (
        'neg_examples.mgf',
        'neg',
        1886,
        {'Cer(t36:0)': (
            MS2Identity(
                score = 51,
                hg = Headgroup(main='Cer', sub=()),
                chainsum = ChainSummary(
                    c = 36,
                    u = 0,
                    typ = ('Sph', 'FA'),
                    attr = (
                        ChainAttr(sph='t', ether=False, oh=()),
                        ChainAttr(sph='', ether=False, oh=())
                    ),
                ),
                chains = (
                    Chain(
                        c = 21,
                        u = 0,
                        typ = 'Sph',
                        attr = ChainAttr(sph='t', ether=False, oh=()),
                        iso = ()
                    ),
                    Chain(
                        c = 15,
                        u = 0,
                        typ = 'FA',
                        attr = ChainAttr(sph='', ether=False, oh=()),
                        iso = ()
                    )
                ),
            ),
        )}
    ),
    (
        'neg_examples.mgf',
        'neg',
        1841,
        {'Cer(t36:0)': (
            MS2Identity(
                score = 48,
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
                        c = 21,
                        u = 0,
                        typ = 'Sph',
                        attr = ChainAttr(sph='t', ether=False, oh=()),
                        iso = ()
                    ),
                    Chain(
                        c = 15,
                        u = 0,
                        typ = 'FA',
                        attr = ChainAttr(sph='', ether=False, oh=()),
                        iso = ()
                    )
                ),
            ),
        )}
    ),
    (
        'neg_examples.mgf',
        'neg',
        1573,
        {'Cer(d36:1-2OH)': (
            MS2Identity(
                score = 93,
                hg = Headgroup(main='Cer', sub=()),
                chainsum = ChainSummary(
                    c = 36,
                    u = 1,
                    typ = ('Sph', 'FAOH'),
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
                        typ = 'FAOH',
                        attr = ChainAttr(sph='', ether=False, oh=('2OH',)),
                        iso = ()
                    )
                ),
            ),
        ),
        'Cer(t36:1)': (
            False,
            MS2Identity(
                score = 38,
                hg = Headgroup(main='Cer', sub=()),
                chainsum = ChainSummary(
                    c = 36,
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
                        u = 1,
                        typ = 'Sph',
                        attr = ChainAttr(sph='t', ether=False, oh=()),
                        iso = ()
                    ),
                    Chain(
                        c = 18,
                        u = 0,
                        typ = 'FAOH',
                        attr = ChainAttr(sph='', ether=False, oh=()),
                        iso = ()
                    )
                ),
            ),
        )}
    ),
    (
        'neg_examples.mgf',
        'neg',
        1596,
        {'Cer(d36:1-2OH)': (
            MS2Identity(
                score = 90,
                hg = Headgroup(main='Cer', sub=()),
                chainsum = ChainSummary(
                    c = 36,
                    u = 1,
                    typ = ('Sph', 'FAOH'),
                    attr = (
                        ChainAttr(sph='d', ether=False, oh=()),
                        ChainAttr(sph='', ether=False, oh=('2OH',))
                    ),
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
                        typ = 'FAOH',
                        attr = ChainAttr(sph='', ether=False, oh=('2OH',)),
                        iso = ()
                    )
                ),
            ),
        ),
        'Cer(t36:1)': (
            False,
            MS2Identity(
                score = 15,
                hg = Headgroup(main='Cer', sub=()),
                chainsum = ChainSummary(
                    c = 36,
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
                        u = 1,
                        typ = 'Sph',
                        attr = ChainAttr(sph='t', ether=False, oh=()),
                        iso = ()
                    ),
                    Chain(
                        c = 18,
                        u = 0,
                        typ = 'FAOH',
                        attr = ChainAttr(sph='', ether=False, oh=()),
                        iso = ()
                    )
                ),
            ),
        )}
    ),
    (
        'neg_examples.mgf',
        'neg',
        2174,
        {'Cer-1P(d34:1)': (
            MS2Identity(
                score = 20,
                hg = Headgroup(main='Cer', sub=('1P',)),
                chainsum = ChainSummary(
                    c = 34,
                    u = 1,
                    typ = ('Sph', 'FA'),
                    attr = (
                        ChainAttr(sph='d', ether=False, oh=()),
                        ChainAttr(sph='', ether=False, oh=())
                    ),
                ),
                chains = None,
            ),
        )}
    ),
    (
        'neg_examples.mgf',
        'neg',
        2265,
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
                        ),
                    iso = None
                    ),
                chains = None,
            ),
        )}
    ),
    (
        'neg_examples.mgf',
        'neg',
        3695,
        {'Hex-Cer(DH42:2)': (
            MS2Identity(
                score = 45,
                hg = Headgroup(main='Cer', sub=('Hex',)),
                chainsum = ChainSummary(
                    c = 42,
                    u = 2,
                    typ = ('Sph', 'FA'),
                    attr = (
                        ChainAttr(sph='DH', ether=False, oh=()),
                        ChainAttr(sph='', ether=False, oh=())
                    ),
                ),
                chains = (
                    Chain(
                        c = 18,
                        u = 1,
                        typ = 'Sph',
                        attr = ChainAttr(sph='DH', ether=False, oh=()),
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
            ),
        )}
    ),
    (
        'neg_examples.mgf',
        'neg',
        2468,
        {'Hex-Cer(t34:1)': (
            MS2Identity(
            score = 78,
            hg = Headgroup(main='Cer', sub=('Hex',)),
            chainsum = ChainSummary(
                c = 34,
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
                    c = 16,
                    u = 1,
                    typ = 'FA',
                    attr = ChainAttr(sph='', ether=False, oh=()),
                    iso = ()
                )
            ),
            ),
        )}
    ),
    # PA positive mode [M+NH4]+
    (
        'pos_examples.mgf',
        'pos',
        2670,
        {'PA(36:2)': (
            MS2Identity(
                score = 20,
                hg = Headgroup(main='PA', sub=()),
                chainsum = ChainSummary(
                    c = 36,
                    u = 2,
                    typ = ('FA', 'FA'),
                    attr = (
                        ChainAttr(sph='', ether=False, oh=()),
                        ChainAttr(sph='', ether=False, oh=())
                    ),
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
            ),
        )}
    ),
    (
        'neg_examples.mgf',
        'neg',
        2646,
        {'PA(36:2)': (
            MS2Identity(
                score = 25,
                hg = Headgroup(main='PA', sub=()),
                chainsum = ChainSummary(
                    c = 36,
                    u = 2,
                    typ = ('FA', 'FA'),
                    attr = (
                        ChainAttr(sph='', ether=False, oh=()),
                        ChainAttr(sph='', ether=False, oh=())
                    ),
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
            ),
        )}
    ),
    # tSph negative [M+HCOO]-
    (
        'neg_examples.mgf',
        'neg',
        506,
        {'Sph-1P(d18:1)': (
            MS2Identity(
                score = 50,
                hg = Headgroup(main='Sph', sub=('1P',)),
                chainsum = ChainSummary(
                    c = 18,
                    u = 1,
                    typ = ('Sph',),
                    attr = (ChainAttr(sph='d', ether=False, oh=()),),
                    iso = None
                    ),
                chains = None,
            ),
        )}
    ),
    (
        'neg_examples.mgf',
        'neg',
        327,
        {'Sph-1P(DH18:0)': (
            MS2Identity(
                score = 40,
                hg = Headgroup(main='Sph', sub=('1P',)),
                chainsum = ChainSummary(
                    c = 18,
                    u = 0,
                    typ = ('Sph',),
                    attr = (ChainAttr(sph='DH', ether=False, oh=()),),
                    iso = None
                    ),
                chains = None,
            ),
        )}
    ),
    # DHSM from GLTPD1 + in vitro
    (
        'pos_examples.mgf',
        'pos',
        2006,
        {'SM(DH34:0)': (
            MS2Identity(
                score = 30,
                hg = Headgroup(main='SM', sub=()),
                chainsum = ChainSummary(
                    c = 34,
                    u = 0,
                    typ = ('Sph', 'FA'),
                    attr = (
                        ChainAttr(sph='DH', ether=False, oh=()),
                        ChainAttr(sph='', ether=False, oh=())
                    ),
                    iso = None
                ),
                chains = None,
            ),
        )}
    ),
    (
        'pos_examples.mgf',
        'pos',
        2032,
        {'SM(DH34:0)': (
            MS2Identity(
                score = 27,
                hg = Headgroup(main='SM', sub=()),
                chainsum = ChainSummary(
                    c = 34,
                    u = 0,
                    typ = ('Sph', 'FA'),
                    attr = (
                        ChainAttr(sph='DH', ether=False, oh=()),
                        ChainAttr(sph='', ether=False, oh=())
                    ),
                    iso = None
                ),
                chains = None,
            ),
        )}
    ),
    (
        'pos_examples.mgf',
        'pos',
        1981,
        {'SM(DH34:0)': (
            MS2Identity(
                score = 30,
                hg = Headgroup(main='SM', sub=()),
                chainsum = ChainSummary(
                    c = 34,
                    u = 0,
                    typ = ('Sph', 'FA'),
                    attr = (
                        ChainAttr(sph='DH', ether=False, oh=()),
                        ChainAttr(sph='', ether=False, oh=())
                    ),
                    iso = None
                ),
                chains = None,
            ),
        )}
    ),
    # PG + in vitro FABP1
    (
        'pos_examples.mgf',
        'pos',
        2040,
        {'PG(36:2)': (
            MS2Identity(
                score = 5,
                max_score = 5,
                score_pct = 100,
                hg = Headgroup(main='PG', sub=()),
                chainsum = ChainSummary(
                    c = 36,
                    u = 2,
                    typ = ('FA', 'FA'),
                    attr = (
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
                        c = 18,
                        u = 1,
                        typ = 'FA',
                        attr = ChainAttr(sph='', ether=False, oh=()),
                        iso = ()
                    )
                ),
            ),
        )}
    ),
        (
        'pos_examples.mgf',
        'pos',
        1980,
        {'PG(36:2)': (
            MS2Identity(
                score = 5,
                max_score = 5,
                score_pct = 100,
                hg = Headgroup(main='PG', sub=()),
                chainsum = ChainSummary(
                    c = 36,
                    u = 2,
                    typ = ('FA', 'FA'),
                    attr = (
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
                        c = 18,
                        u = 1,
                        typ = 'FA',
                        attr = ChainAttr(sph='', ether=False, oh=()),
                        iso = ()
                    )
                ),
            ),
        )}
    ),
    # SM + from in vitro GLTPD1
    (
        'pos_examples.mgf',
        'pos',
        1940,
        {'SM(d35:1)': (
            MS2Identity(
                score = 27,
                max_score = 49,
                score_pct = 55,
                hg = Headgroup(main='SM', sub=()),
                chainsum = ChainSummary(
                    c = 35,
                    u = 1,
                    typ = ('Sph', 'FA'),
                    attr = (
                        ChainAttr(sph='d', ether=False, oh=()),
                        ChainAttr(sph='', ether=False, oh=())
                    ),
                    iso = None
                ),
                chains = None,
            ),
        )}
    ),
    # LysoPG + from in vitro SCP2D1
    (
        'pos_examples.mgf',
        'pos',
        553,
        {'Lyso-PG(18:1)': (
            MS2Identity(
                score = 5,
                max_score = 10,
                score_pct = 50,
                hg = Headgroup(main='PG', sub=('Lyso',)),
                chainsum = ChainSummary(
                    c = 18,
                    u = 1,
                    typ = ('FA',),
                    attr = (ChainAttr(sph='', ether=False, oh=()),),
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
                ),
            ),
        )}
    ),
    # LysoPC + from in vitro SCP2D1
    (
        'pos_examples.mgf',
        'pos',
        436,
        {'Lyso-PC(20:4)': (
            MS2Identity(
                score = 10,
                max_score = 15,
                score_pct = 67,
                hg = Headgroup(main='PC', sub=('Lyso',)),
                chainsum = ChainSummary(
                    c = 20,
                    u = 4,
                    typ = ('FA',),
                    attr = (ChainAttr(sph='', ether=False, oh=()),),
                    iso = None
                ),
                chains = None,
            ),
        )}
    ),
    # FA(18:1) negative from SCP2D1 in vitro
    (
        'neg_examples.mgf',
        'neg',
        826,
        {'FA(18:1)': (
            MS2Identity(
                score = 10,
                max_score = 10,
                score_pct = 100,
                hg = Headgroup(main='FA', sub=()),
                chainsum = ChainSummary(
                    c = 18,
                    u = 1,
                    typ = ('FA',),
                    attr = (ChainAttr(sph='', ether=False, oh=()),),
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
                ),
            ),
        )}
    ),
    # LysoPE negative from SCP2D1 in vitro
    (
        'neg_examples.mgf',
        'neg',
        343,
        {'Lyso-PE(16:1)': (
            MS2Identity(
                score = 8,
                max_score = 11,
                score_pct = 73,
                hg = Headgroup(main='PE', sub=('Lyso',)),
                chainsum = ChainSummary(
                    c = 16,
                    u = 1,
                    typ = ('FA',),
                    attr = (ChainAttr(sph='', ether=False, oh=()),),
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
                ),
            ),
        )}
    ),
    # PS negative from ORP9 in vitro
    (
        'neg_examples.mgf',
        'neg',
        1710,
        {'PS(36:1)': (
            MS2Identity(
                score = 11,
                max_score = 20,
                score_pct = 55,
                hg = Headgroup(main='PS', sub=()),
                chainsum = ChainSummary(
                    c = 36,
                    u = 1,
                    typ = ('FA', 'FA'),
                    attr = (
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
                        c = 18,
                        u = 0,
                        typ = 'FA',
                        attr = ChainAttr(sph='', ether=False, oh=()),
                        iso = ()
                    )
                ),
            ),
        )}
    ),
    # PS positive from BPI in vivo
    (
        'pos_examples.mgf',
        'pos',
        2403,
        {'PS(34:1)': (
            MS2Identity(
                score = 5,
                max_score = 5,
                score_pct = 100,
                hg = Headgroup(main='PS', sub=()),
                chainsum = ChainSummary(
                    c = 34,
                    u = 1,
                    typ = ('FA', 'FA'),
                    attr = (
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
                        c = 16,
                        u = 0,
                        typ = 'FA',
                        attr = ChainAttr(sph='', ether=False, oh=()),
                        iso = ()
                    )
                ),
            ),
        )}
    ),
    # PI positive BPI in vivo
    (
        'pos_examples.mgf',
        'pos',
        2397,
        {'PI(34:1)': (
            MS2Identity(
                score = 5,
                max_score = 9,
                score_pct = 56,
                hg = Headgroup(main='PI', sub=()),
                chainsum = ChainSummary(
                    c = 34,
                    u = 1,
                    typ = ('FA', 'FA'),
                    attr = (
                        ChainAttr(sph='', ether=False, oh=()),
                        ChainAttr(sph='', ether=False, oh=())
                    ),
                    iso = None
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
            ),
        )}
    ),
]


class TestScan(object):
    """ """
    
    @pytest.mark.parametrize('mgfname, ionmode, scan_id, identity', specimens)
    def test_identification(
            self,
            mgfname,
            ionmode,
            scan_id,
            identity
        ):
        """

        Parameters
        ----------
        mgfname :
            
        ionmode :
            
        scan_id :
            
        identity :
            

        Returns
        -------

        """
        
        std_key = (scan_id, ionmode)
        
        scan_args = (
            special_settings[std_key]['scan_args']
            if (
                std_key in special_settings and
                'scan_args' in special_settings[std_key]
            )
            else {}
        )
        
        mgfpath = os.path.join(common.ROOT, 'data', 'ms2_examples', mgfname)
        scan = ms2.Scan.from_mgf(mgfpath, scan_id, ionmode, **scan_args)
        
        iscan = scan.identify()
        
        sm_rec = list(
            scan.get_ms1_records(hg = 'SM', databases = {'lipyd.lipid'})
        )
        
        highest_score = max(
            (
                i.score_pct
                for ids in iscan.values()
                for i in ids
            ),
            default = 0,
        )
        
        for name, ids in identity.items():
            
            assert name in iscan
            
            false_identity = ids[0] == False
            
            if false_identity:
                
                ids = ids[1:]
            
            for i in ids:
                
                assert i in iscan[name]
            
            highest_for_name = max(
                (i.score_pct for i in iscan[name]),
                default = 0,
            )
            
            assert (
                highest_for_name == highest_score or (
                    false_identity and
                    highest_for_name < highest_score
                )
            )
