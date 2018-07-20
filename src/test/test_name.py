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

import lipyd.name
import lipyd.lipproc


swl_names = {
    'Phosphatidylcholine(36:1)':
        (
            lipyd.lipproc.Headgroup(main = 'PC', sub = ('',)),
            lipyd.lipproc.Chain(
                c = 36, u = 1, typ = None, iso = (),
                attr = lipyd.lipproc.ChainAttr(
                    sph = '', ether = False, oh = ()
                )
            ),
            []
        ),
    'Phosphatidylinositol(18:1/18:0)':
        (
            lipyd.lipproc.Headgroup(main = 'PI', sub = ('',)),
            lipyd.lipproc.Chain(
                c = 36, u = 1, typ = None, iso = (),
                attr = lipyd.lipproc.ChainAttr(
                    sph = '', ether = False, oh = ()
                )
            ),
            (
                lipyd.lipproc.Chain(
                    c = 18, u = 1, typ = 'FA', iso = (),
                    attr = lipyd.lipproc.ChainAttr(
                        sph = '', ether = False, oh = ()
                    )
                ),
                lipyd.lipproc.Chain(
                    c = 18, u = 0, typ = 'FA', iso = (),
                    attr = lipyd.lipproc.ChainAttr(
                        sph = '', ether = False, oh = ()
                    )
                )
            )
        )
}

lmp_names = {
    'Cer(t18:1/18:1)':
        (
            lipyd.lipproc.Headgroup(main='Cer', sub=('',)),
            lipyd.lipproc.Chain(
                c = 36, u = 2, typ = None, iso = (),
                attr = lipyd.lipproc.ChainAttr(
                    sph = 't', ether = False, oh = ()
                )
            ),
            (
                lipyd.lipproc.Chain(
                    c = 18, u = 1, typ = 'Sph', iso = (),
                    attr = lipyd.lipproc.ChainAttr(
                        sph = 't', ether = False, oh = ()
                    )
                ),
                lipyd.lipproc.Chain(
                    c = 18, u = 1, typ = 'FA', iso = (),
                    attr = lipyd.lipproc.ChainAttr(
                        sph = '', ether = False, oh = ()
                    )
                )
            )
        ),
    'Cer(d18:0/18:1)':
        (
            lipyd.lipproc.Headgroup(main = 'Cer', sub = ('',)),
            lipyd.lipproc.Chain(
                c = 36, u = 1, typ = None, iso = (),
                attr = lipyd.lipproc.ChainAttr(
                    sph = 'DH', ether = False, oh = ()
                )
            ),
            (
                lipyd.lipproc.Chain(
                    c = 18, u = 0, typ = 'Sph', iso = (),
                    attr = lipyd.lipproc.ChainAttr(
                        sph = 'DH', ether = False, oh = ()
                    )
                ),
                lipyd.lipproc.Chain(
                    c = 18, u = 1, typ = 'FA', iso = (),
                    attr = lipyd.lipproc.ChainAttr(
                        sph = '', ether = False, oh = ())
                )
            )
        )
}


class TestName(object):
    
    nameproc = lipyd.name.LipidNameProcessor()
    
    def test_name_swisslipids(self):
        
        for name, result in swl_names.items():
            
            assert (
                self.nameproc.process(name, database = 'swisslipids') ==
                result
            )
    
    def test_name_lipidmaps(self):
        
        for name, result in lmp_names.items():
            
            assert (
                self.nameproc.process(name, database = 'lipidmaps') ==
                result
            )
    
    def test_name_iso(self):
        
        self.nameproc.iso = True
        self.nameproc.database = 'lipidmaps'
        
        pelghl = self.nameproc.process(
            'PE(16:0/18:1(9Z))-15-isoLG hydroxylactam'
        )
        
        pelghl_hg = lipyd.lipproc.Headgroup(main = 'PE', sub = ('LGHL',))
        pelghl_chainsum = lipyd.lipproc.Chain(
            c = 34, u = 1, typ = None, iso = (),
            attr = lipyd.lipproc.ChainAttr(sph = '', ether = False, oh = ())
        )
        pelghl_chains = (
            lipyd.lipproc.Chain(
                c = 16, u = 0, typ = 'FA', iso = (),
                attr = lipyd.lipproc.ChainAttr(
                    sph = '', ether = False, oh = ()
                )
            ),
            lipyd.lipproc.Chain(
                c = 18, u = 1, typ = 'FA', iso = ('9Z',),
                attr = lipyd.lipproc.ChainAttr(
                    sph = '', ether = False, oh = ()
                )
            )
        )
        
        assert pelghl_hg == pelghl[0]
        assert pelghl_chainsum == pelghl[1]
        assert pelghl_chains == pelghl[2]
        
        fahfa = self.nameproc.process('FAHFA(16:0/10-O-18:0)')
        
        fahfa_hg = lipyd.lipproc.Headgroup(main='FAHFA', sub=('',))
        fahfa_chainsum = lipyd.lipproc.Chain(
            c = 34, u = 0, typ = None, iso = (),
            attr = lipyd.lipproc.ChainAttr(sph='', ether=False, oh=())
        )
        fahfa_chains = (
            lipyd.lipproc.Chain(
                c = 16, u = 0, typ = 'FA', iso = (),
                attr = lipyd.lipproc.ChainAttr(
                    sph = '', ether = False, oh = ()
                )
            ),
            lipyd.lipproc.Chain(
                c = 18, u = 0, typ = 'FA', iso = (),
                attr = lipyd.lipproc.ChainAttr(
                    sph = '', ether = False, oh = ()
                )
            )
        )
        
        assert fahfa_hg == fahfa[0]
        assert fahfa_chainsum == fahfa[1]
        assert fahfa_chains == fahfa[2]
