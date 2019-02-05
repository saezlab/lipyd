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
#  Website: http://denes.omnipathdb.org/
#

import re
import collections
import itertools
import functools


FAMILIES = {
    'SL':  {'Sph', 'Cer', 'SM'},
    'GL':  {'MAG', 'DAG', 'TAG'},
    'GPL': {
        'PC', 'PE', 'PE', 'PG', 'BMP',
        'PA', 'PI', 'PIP', 'PIP2', 'PIP3'
    },
    'GLGPL': {
        'MAG', 'DAG', 'TAG',
        'PC', 'PE', 'PE', 'PG', 'BMP',
        'PA', 'PI', 'PIP', 'PIP2', 'PIP3'
    }
}


SUBCLS_PRE = {'Hex', 'Hex2', 'SHex', 'SHex2', 'Lyso'}

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
    """ """
    
    def __new__(cls, c, u, typ = 'FA', attr = ChainAttr(), iso = ()):
        
        if hasattr(attr, 'sph') and attr.sph == 'd' and u == 0:
            
            attr = ChainAttr(sph = 'DH', ether = attr.ether, oh = attr.oh)
        
        return super(Chain, cls).__new__(
            cls, c, u, typ = typ, attr = attr, iso = iso
        )
    
    def __str__(self, iso = False):
        
        return '%s%s%u:%u%s%s' % (
            # ether prefix
            'O-' if self.attr.ether else '',
            # sphingoid base prefix e.g. dCer(d38:1)
            self.attr.sph,
            # the carbon count
            self.c,
            # the unsaturation
            self.u,
            # isomer information
            '(%s)' % ','.join(self.iso) if self.iso else '',
            # postfix for hydroxylated fatty acyls e.g. PC(32:1-2OH)
            '-%s' % '-'.join(self.attr.oh) if self.attr.oh else ''
        )
    
    def isomer_str(self):
        """ """
        
        return self.__str__(iso = True)
    
    def __add__(self, other):
        
        return sum_chains((self, other))


class ChainSummary(Chain):
    """ """
    
    def __new__(cls, c, u, typ = (), attr = (), iso = None):
        
        if u == 0 and attr and attr[0].sph == 'd':
            
            attr = (
                (
                    # fisrt chain, sphingosine base
                    ChainAttr(
                        sph = 'DH',
                        ether = attr[0].ether,
                        oh = attr[0].oh
                    ),
                ) +
                # all other chains
                attr[1:]
            )
        
        return super(ChainSummary, cls).__new__(
            cls, c, u, typ = typ, attr = attr, iso = None
        )
    
    def __str__(self):
        
        return '%s%s%u:%u%s' % (
            # ether prefix
            'O-' if any(a.ether for a in self.attr) else '',
            # sphingoid base prefix e.g. dCer(d38:1)
            ''.join(a.sph for a in self.attr), # of course max one of these
                                               # is not empty string
            # the carbon count
            self.c,
            # the unsaturation
            self.u,
            # postfix for hydroxylated fatty acyls e.g. PC(32:1-2OH)
            '-%s' % '-'.join('-'.join(a.oh) for a in self.attr if a.oh)
                if any(a.oh for a in self.attr) else ''
        )
    
    def __len__(self):
        
        return len(self.typ)


Headgroup = collections.namedtuple(
    'Headgroup',
    ['main', 'sub']
)
Headgroup.__new__.__defaults__ = ((),)


LipidLabel = collections.namedtuple(
    'LipidLabel',
    ['db_id', 'db', 'names', 'formula']
)
# names are empty tuple by default
LipidLabel.__new__.defaults__ = ((), None)


class LipidRecord(collections.namedtuple(
        'LipidRecordBase',
        ['lab', 'hg', 'chainsum', 'chains']
    )):
    """ """
    
    def __new__(cls, lab, hg, chainsum, chains):
        
        return super(LipidRecord, cls).__new__(cls, lab, hg, chainsum, chains)
    
    def full_str(self):
        """ """
        
        return full_str(self.hg, self.chains, iso = False)
    
    def summary_str(self):
        """ """
        
        return summary_str(self.hg, self.chainsum)
    
    def subclass_str(self):
        """ """
        
        return subclass_str(self.hg, self.chainsum)


def empty_chain():
    """Returns an empty Chain object which might serve as a dummy object."""
    
    return Chain(c = 0, u = 0, attr = ChainAttr())


def empty_chainsum():
    """Returns an empty ChainSummary object."""
    
    return ChainSummary(c = 0, u = 0, attr = (), typ = ())


def str2hg(hgstr):
    """From a headgroup string representation creates a Headgroup object.

    Parameters
    ----------
    hgstr :
        

    Returns
    -------

    """
    
    pieces = hgstr.split('-')
    
    return Headgroup(main = pieces[-1], sub = tuple(pieces[:-1]))


def sum_chains(chains):
    """From a list of chains creates a summary Chain object.

    Parameters
    ----------
    chains :
        

    Returns
    -------

    """
    
    return empty_chainsum() if not chains else (
        ChainSummary(
            c = sum(i.c for i in chains),
            u = sum(i.u for i in chains),
            attr = tuple(c.attr for c in chains),
            typ = tuple(c.typ for c in chains)
        )
    )


def collapse_attrs(chains):
    """Combine the attributes of arbitrary number of chains.

    Parameters
    ----------
    chains :
        

    Returns
    -------

    """
    
    return functools.reduce(combine_attrs, (c.attr for c in chains))


def combine_attrs(a1, a2):
    """Combines the attributes of 2 chains.

    Parameters
    ----------
    a1 :
        
    a2 :
        

    Returns
    -------

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

    Parameters
    ----------
    hg : lipyd.lipproc.Headgroup
        ``Headgroup`` object.
    chainsum : lipyd.lipproc.ChainSummary
        ``ChainSummary`` object
    
    Returns
    -------
    String representation with sum of chain lengths and unsaturations.
    """
    
    subcls_pre, sphingo_prefix, ether_prefix, subcls_post, hydroxy = (
        get_attributes(hg, chainsum)
    )
    
    return '%s%s%s%s' % (
        # subclass attributes like *PE*-Cer, *Lyso*-PC
        subcls_pre,
        # main class of headgroup e.g. Cer, PS
        hg.main,
        # subclass attributes like 1-O-phosphate group of Cer1P, Sph1P, etc
        subcls_post,
        # chains summary
        ('(%s)' % chainsum.__str__())
            if chainsum is not None and chainsum.c > 0 else
        ''
    )


def full_str(hg, chains, iso = False):
    """
    From a Headgroup and a tuple of Chain objects returns a

    Parameters
    ----------
    hg : Headgroup
        ``Headgroup`` object.
    chains : tuple
        Tuple of ``Chain`` objects.
    iso :
        Include isomer information.

    Returns
    -------
    String representation with information about chain lengths and
    unsaturation, optionally E/Z isomerism.
    """
    
    subcls_pre, sphingo_prefix, ether_prefix, subcls_post, hydroxy = (
        get_attributes(hg, sum_chains(chains))
    )
    
    return '%s%s%s%s' % (
        # subclass attributes like *PE*-Cer, *Lyso*-PC
        subcls_pre,
        # main class of headgroup e.g. Cer, PS
        hg.main,
        # subclass attributes like 1-O-phosphate group of Cer1P, Sph1P, etc
        subcls_post,
        # chains
        ('(%s)' % '/'.join(c.__str__(iso = iso) for c in chains))
            if chains else
        ''
    )

def subclass_str(hg, chainsum = None):
    """
    From Headgroup and summary Chain object creates a subclass level
    headgroup string.

    Parameters
    ----------
    hg :
        
    chainsum :
         (Default value = None)

    Returns
    -------

    """
    
    subcls_pre, sphingo_prefix, ether_prefix, subcls_post, hydroxy = (
        get_attributes(hg, chainsum)
    )
    
    return '%s%s%s%s%s%s' % (
        # subclass attributes like *PE*-Cer, *Lyso*-PC
        subcls_pre,
        # prefix of shingoid base subclass: d, t, k, DH
        sphingo_prefix,
        # main class of headgroup e.g. Cer, PS
        hg.main,
        # postfix of ether lipids e.g. PC-O, PE-O
        '-O' if ether_prefix else '',
        # subclass attributes like 1-O-phosphate group of Cer1P, Sph1P, etc
        subcls_post,
        # postfix of hydroxylated fatty acyl e.g. Cer-2OH
        hydroxy
    )

def class_str(hg):
    """From Headgroup returns the main class as string.

    Parameters
    ----------
    hg :
        

    Returns
    -------

    """
    
    return hg.main

def get_attributes(hg, chainsum = None):
    """Processes a Headgroup and a summary Chain object and returns the
    name pre- and postfix string elements.

    Parameters
    ----------
    hg :
        
    chainsum :
         (Default value = None)

    Returns
    -------

    """
    
    chainsum = chainsum or empty_chainsum()
    
    hydroxy = '-'.join('-'.join(c.oh) for c in chainsum.attr)
    hydroxy = '-%s' % hydroxy if hydroxy else ''
    
    subcls_pre   = '-'.join(i for i in hg.sub if i in SUBCLS_PRE)
    subcls_pre   = '%s-' % subcls_pre if subcls_pre else ''
    subcls_post  = '-'.join(i for i in hg.sub if i not in SUBCLS_PRE)
    subcls_post  = '-%s' % subcls_post if subcls_post else ''
    
    sphingo_prefix = ''.join(a.sph for a in chainsum.attr)
    ether_prefix = any(a.ether for a in chainsum.attr)
    # this I leave here if maybe later I decide to have
    # plain tuples...
    #
    # sphingo_prefix = set(chainsum.p) & {'d', 't', 'k', 'DH'}
    # sphingo_prefix = sphingo_prefix.pop() if sphingo_prefix else ''
    # ether_prefix = 'O' in chainsum.p
    
    return subcls_pre, sphingo_prefix, ether_prefix, subcls_post, hydroxy


def match_constraint(rec, constr):
    """Matches an MS2 fragment constraint (fragment.FragConstraint)
    against a lipid record lipproc.LipidRecord.
    
    Returns the indices of the chains as integers in a set.
    
    The default attribute values of the `FragConstraint` object always
    correspond to bypass all filters here. Defining a specific value
    may limit the records complying with the constraint.
    
    Args
    ----

    Parameters
    ----------
    LipidRecord :
        rec:
        A lipid database record object.
    fragment :
        FragConstraint constr:
        An MS2 fragment constraint object.
    rec :
        
    constr :
        

    Returns
    -------

    """
    
    match = False
    chains = set()
    
    if (
        (
            constr.hg == rec.hg.main or (
                constr.hg is None and
                constr.family in FAMILIES and
                rec.hg.main in FAMILIES[constr.family]
            ) or (
                constr.hg is None and
                constr.family is None
            )
        ) and (
            constr.sub is None or
            set(constr.sub) == set(rec.hg.sub)
        )
    ):
        
        match = True
        chainsum = rec.chainsum if rec.chainsum else sum_chains(rec.chains)
        
        for (i, attr), rec_chaintype in zip(
            enumerate(chainsum.attr), chainsum.typ
        ):
            
            if (
                (
                    constr.chaintype is None or
                    constr.chaintype == rec_chaintype
                ) and (
                    constr.sph is None or
                    constr.sph == attr.sph
                ) and
                # matching only the number of OH groups
                (
                    len(attr.oh) == constr.oh
                    if type(constr.oh) is int else
                    set(attr.oh) == set(constr.oh)
                )
            ):
                
                chains.add(i)
    
    return match, chains


def match_constraints(rec, constraints):
    """Matches all fragment constraints in the iterable `constraints`
    against all chains in MS1 record `rec`.
    
    Returns a boolean (wether the fragment can possibly origin from the
    molecular species in the record) and a tuple of chain positions which
    can be the source of the fragment if the fragment is from an aliphatic
    chain moiety.
    
    Args
    ----

    Parameters
    ----------
    LipidRecord :
        rec:
        An MS1 database record object.
    iterable :
        constraints:
        A number of `fragment.FragConstraint` objects.
    rec :
        
    constraints :
        

    Returns
    -------

    """
    
    match = False
    chains = set()
    
    for constr in constraints:
        
        match, chains_ = match_constraint(rec, constr)
        chains.update(chains_)
    
    return match, chains

def cu_str(c, u):
    """

    Parameters
    ----------
    c :
        
    u :
        

    Returns
    -------

    """
    
    return '%u:%u' % (c, u)

def charge_str(charge):
    """

    Parameters
    ----------
    charge :
        

    Returns
    -------

    """
    
    return (
        '' if charge == 0 else
        '-' if charge < 0 else
        '+'
    )

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
    r'(?:\(?([0-9EZ,]{2,})\)?)?'
    r'(?:[-\(]([0-9]{0,2}OH)\)?)?'
    r'[/_]?'
    # 2
    r'((?:[0-9]+-)?'
    r'[POdtDH]{0,2})-?'
    r'([0-9]{0,2}):?([0-9]{0,2})'
    r'(?:\(?([0-9EZ,]{2,})\))?'
    r'(?:[-\(]([0-9]{0,2}OH)\)?)?'
    r'[/_]?'
    # 3
    r'((?:[0-9]+-)?'
    r'[POdtDH]{0,2})-?'
    r'([0-9]{0,2}):?([0-9]{0,2})'
    r'(?:\(?([0-9EZ,]{2,})\))?'
    r'(?:[-\(]([0-9]{0,2}OH)\)?)?'
    r'[/_]?'
    # 4
    r'((?:[0-9]+-)?'
    r'[POdtDH]{0,2})-?'
    r'([0-9]{0,2}):?([0-9]{0,2})'
    r'(?:\(?([0-9EZ,]*)\)?)?'
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
