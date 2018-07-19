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

import re
import collections
import itertools
import functools


ChainAttr = collections.namedtuple(
    'ChainAttr',
    ['sph', 'ether', 'oh']
)
# attributes of a chain e.g. d or t Ceramide, 2OH fatty acyl, etc
# defaults:
ChainAttr.__new__.__defaults__ = ('', False, ())


class Chain(collections.namedtuple(
        'ChainBase',
        ['c', 'u', 'typ', 'attr', 'iso']
    )):
    
    def __new__(cls, c, u, typ = 'FA', attr = ChainAttr(), iso = ()):
        
        if attr.sph == 'd' and u == 0:
            
            attr = ChainAttr(sph = 'DH', ether = attr.ether, oh = attr.oh)
        
        return super(Chain, cls).__new__(
            cls, c, u, typ = typ, attr = attr, iso = iso
        )


Headgroup = collections.namedtuple(
    'Headgroup',
    ['main', 'sub']
)
Headgroup.__new__.__defaults__ = ((),)


LipidLabel = collections.namedtuple(
    'LipidLabel',
    ['cls', 'db_id', 'db', 'names']
)
# names are empty tuple by default
LipidLabel.__new__.defaults__ = ((),)


def empty_chain():
    """
    Returns an empty Chain object which might serve as a dummy object.
    """
    
    return Chain(c = 0, u = 0, attr = ChainAttr())


def str2hg(hgstr):
    """
    From a headgroup string representation creates a Headgroup object.
    """
    
    pieces = hgstr.split('-')
    
    return Headgroup(main = pieces[-1], sub = tuple(pieces[:-1]))


def sum_chains(chains):
    """
    From a list of chains creates a summary Chain object.
    """
    
    return (
        Chain(
            c = sum(i.c for i in chains),
            u = sum(i.u for i in chains),
            attr = collapse_attrs(chains),
            typ = None
        )
    )


def collapse_attrs(chains):
    """
    Combine the attributes of arbitrary number of chains.
    """
    
    return functools.reduce(combine_attrs, (c.attr for c in chains))


def combine_attrs(a1, a2):
    """
    Combines the attributes of 2 chains.
    """
    
    return ChainAttr(
        sph = a1.sph or a2.sph,
        ether = a1.ether or a2.ether,
        oh = tuple(itertools.chain(a1.oh, a2.oh))
    )

def summary_str(hg, chainsum):
    """
    Creates a summary string representation from the headgroup name and
    a summary Chain object.
    """
    
    subcls, sphingo_prefix, ether_prefix, p1, hydroxy = get_attributes(
        hg, chaisum
    )
    
    return '%s%s%s(%s%s%u:%u%s)' % (
        # subclass attributes like *PE*-Cer, *Lyso*-PC
        subcls,
        # main class of headgroup e.g. Cer, PS
        hg.main,
        # 1-O-phosphate group of Cer1P, Sph1P, etc
        p1,
        # prefix of carbon counts at ether lipids e.g. PC(O-36:1)
        'O-' if ether_prefix else '',
        # sphingoid base prefix e.g. dCer(d38:1)
        sphingo_prefix,
        # the total carbon count
        chainsum.c,
        # the total unsaturation
        chainsum.u,
        # postfix for hydroxylated fatty acyls e.g. PC(32:1-2OH)
        hydroxy
    )


def subclass(hg, chainsum = None):
    """
    From Headgroup and summary Chain object creates a subclass level
    headgroup string.
    """
    
    subcls, sphingo_prefix, ether_prefix, p1, hydroxy = get_attributes(
        hg, chaisum
    )
    
    return '%s%s%s%s%s%s' % (
        # subclass attributes like *PE*-Cer, *Lyso*-PC
        subcls,
        # prefix of shingoid base subclass: d, t, k, DH
        sphingo_prefix,
        # main class of headgroup e.g. Cer, PS
        hg.main,
        # postfix of ether lipids e.g. PC-O, PE-O
        '-O' if ether_prefix else '',
        # 1-O-phosphate group of Cer1P, Sph1P, etc
        p1,
        # postfix of hydroxylated fatty acyl e.g. Cer-2OH
        hydroxy
    )


def get_attributes(hg, chaisum = None):
    """
    Processes a Headgroup and a summary Chain object and returns the
    name pre- and postfix string elements.
    """
    
    chainsum = chainsum or empty_chain()
    
    hydroxy = '-'.join(chainsum.attr.hydroxy)
    hydroxy = '-%s' % hydroxy if hydroxy else ''
    
    subcls  = '-'.join(i for i in hg.sub if i != '1P')
    subcls  = '%s-' % subcls if subcls else ''
    p1 = '1P' if '1P' in hg.sub else ''
    
    sphingo_prefix = chainsum.attr.sph
    ether_prefix = chainsum.attr.ether
    # this I leave here if maybe later I decide to have
    # plain tuples...
    #
    # sphingo_prefix = set(chainsum.p) & {'d', 't', 'k', 'DH'}
    # sphingo_prefix = sphingo_prefix.pop() if sphingo_prefix else ''
    # ether_prefix = 'O' in chainsum.p
    
    return subcls, sphingo_prefix, ether_prefix, p1, hydroxy


# regex captures the summary carbon count
rechainsum = re.compile(
    r'\('
    # prefix (d, t , DH, O-, P-)
    r'([POdtDH]{0,2})-?'
    # cc and unsat
    r'([0-9]{1,2}):([0-9]{1,2})'
    # optional OH
    r'(?:[-\(]([0-9]{0,2}OH)\)?)?'
    r'\)'
)

# captures 1-4 aliphatic chains data
rechain = re.compile(
    r'\('
    # 1
    r'([POdtDH]{0,2})-?'
    r'([0-9]{1,2}):([0-9]{1,2})'
    r'(?:[-\(]([0-9]{0,2}OH)\)?)?'
    r'[/_]?'
    # 2
    r'([POdtDH]{0,2})-?'
    r'([0-9]{0,2}):?([0-9]{0,2})'
    r'(?:[-\(]([0-9]{0,2}OH)\)?)?'
    r'[/_]?'
    # 3
    r'([POdtDH]{0,2})-?'
    r'([0-9]{0,2}):?([0-9]{0,2})'
    r'(?:[-\(]([0-9]{0,2}OH)\)?)?'
    r'[/_]?'
    # 4
    r'([POdtDH]{0,2})-?'
    r'([0-9]{0,2}):?([0-9]{0,2})'
    r'(?:[-\(]([0-9]{0,2}OH)\)?)?'
    r'\)'
)

# captures 1-4 aliphatic chains with
# conformational isomeric information
rechainiso = re.compile(
    r'\(?'
    # 1
    r'((?:[0-9]+-)?'
    r'[POdtDH]{0,2})-?'
    r'([0-9]{1,2}):([0-9]{1,2})'
    r'\(?([0-9EZ,]*)\)?'
    r'(?:[-\(]([0-9]{0,2}OH)\)?)?'
    r'[/_]?'
    # 2
    r'((?:[0-9]+-)?'
    r'[POdtDH]{0,2})-?'
    r'([0-9]{0,2}):?([0-9]{0,2})'
    r'\(?([0-9EZ,]*)\)?'
    r'(?:[-\(]([0-9]{0,2}OH)\)?)?'
    r'[/_]?'
    # 3
    r'((?:[0-9]+-)?'
    r'[POdtDH]{0,2})-?'
    r'([0-9]{0,2}):?([0-9]{0,2})'
    r'\(?([0-9EZ,]*)\)?'
    r'(?:[-\(]([0-9]{0,2}OH)\)?)?'
    r'[/_]?'
    # 4
    r'((?:[0-9]+-)?'
    r'[POdtDH]{0,2})-?'
    r'([0-9]{0,2}):?([0-9]{0,2})'
    r'\(?([0-9EZ,]*)\)?'
    r'(?:[-\(]([0-9]{0,2}OH)\)?)?'
    r'\)?'
)

# methyl or ethyl
reme = re.compile(r'methyl|ethyl')

# ?
rebr = re.compile(
    r'(1(?:,2-di)?)-\(((?:[0-9]{0,2}[-]?methyl|ethyl)?)[A-z0-9-]+\)'
    r'-([2,3]{1,3}(?:-di)?)-'
    r'\(((?:[0-9]{0,2}[-]?methyl|ethyl)?)[A-z0-9-]+\)'
)
