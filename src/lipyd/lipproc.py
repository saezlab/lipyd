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


Chain = collections.namedtuple(
    'Chain',
    ['c', 'u', 't', 'p', 'i']
)
# by default type id `FA` the prefix is empty tuple
Chain.__new__.__defaults__ = ('FA', (), ())

LipidLabel = collections.namedtuple(
    'LipidLabel',
    ['cls', 'db_id', 'db', 'names']
)
# names are empty tuple by default
LipidLabel.__new__.defaults__ = ((),)

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
    r'([POdtDH]{0,2})-?'
    r'([0-9]{1,2}):([0-9]{1,2})'
    r'(?:[-\(]([0-9]{0,2}OH)\)?)?'
    r'[/_]?'
    r'([POdtDH]{0,2})-?'
    r'([0-9]{0,2}):?([0-9]{0,2})'
    r'(?:[-\(]([0-9]{0,2}OH)\)?)?'
    r'[/_]?'
    r'([POdtDH]{0,2})-?'
    r'([0-9]{0,2}):?([0-9]{0,2})'
    r'(?:[-\(]([0-9]{0,2}OH)\)?)?'
    r'[/_]?'
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
    r'((?:[-\(]2OH\)?)?)'
    r'[/_]?'
    # 2
    r'((?:[0-9]+-)?'
    r'[POdtDH]{0,2})-?'
    r'([0-9]{0,2}):?([0-9]{0,2})'
    r'\(?([0-9EZ,]*)\)?'
    r'((?:[-\(]2OH\)?)?)'
    r'[/_]?'
    # 3
    r'((?:[0-9]+-)?'
    r'[POdtDH]{0,2})-?'
    r'([0-9]{0,2}):?([0-9]{0,2})'
    r'\(?([0-9EZ,]*)\)?'
    r'((?:[-\(]2OH\)?)?)'
    r'[/_]?'
    # 4
    r'((?:[0-9]+-)?'
    r'[POdtDH]{0,2})-?'
    r'([0-9]{0,2}):?([0-9]{0,2})'
    r'\(?([0-9EZ,]*)\)?'
    r'((?:[-\(]2OH\)?)?)'
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
