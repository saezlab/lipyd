#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `emese` python module
#
#  Copyright (c) 2015-2017 - EMBL-EBI
#
#  File author(s): Dénes Türei (turei.denes@gmail.com)
#
#  This code is not for public use.
#  Please do not redistribute.
#  For permission please contact me.
#
#  Website: http://www.ebi.ac.uk/~denes
#

from __future__ import print_function
from past.builtins import xrange, range, reduce

import sys
import warnings
import traceback
import itertools

import emese.settings as settings

if 'unicode' not in __builtins__:
    unicode = str

simpleTypes = [int, float, str, unicode]

charTypes = [str, unicode]

def uniqList(seq):
    # Not order preserving
    # from http://www.peterbe.com/plog/uniqifiers-benchmark
    keys = {}
    for e in seq:
        try:
            keys[e] = 1
        except:
            sys.stdout.write(e)
            sys.stdout.write('\n')
            sys.stdout.write(seq)
            sys.stdout.write('\n')
            sys.stdout.write(keys)
            sys.stdout.write('\n')
    return keys.keys()

def flatList(lst):
    return [it for sl in lst for it in sl]

def delEmpty(lst):
    return [i for i in lst if len(i) > 0]

def uniqOrdList(seq, idfun = None): 
   # Order preserving
   # from http://www.peterbe.com/plog/uniqifiers-benchmark
   if idfun is None:
       def idfun(x): return x
   seen = {}
   result = []
   for item in seq:
       marker = idfun(item)
       if marker in seen: continue
       seen[marker] = 1
       result.append(item)
   return result

def addToList(lst, toadd):
    if isinstance(toadd, list):
        lst += toadd
    else:
        lst.append(toadd)
    if None in lst:
        lst.remove(None)
    return uniqList(lst)

def something(anything):
    return not (anything is None or \
        (type(anything) in [list, set, dict, str, unicode] \
            and len(anything) == 0))

def warn_with_traceback(message, category, filename, lineno, file=None, line=None):
    """
    Prints warnings with traceback.
    From https://stackoverflow.com/a/22376126/854988
    """
    traceback.print_stack()
    log = file if hasattr(file,'write') else sys.stderr
    log.write(warnings.formatwarning(message, category, filename, lineno, line))

def get_param(par):
    
    if hasattr(settings, par):
        return getattr(settings, par)

def set_param(par, val):
    
    setattr(settings, par, val)

fa_greek_parts = {
    'cc': {
        'hex': 6,
        'hept': 7,
        'oct': 8,
        'non': 9,
        'dec': 10,
        'undec': 11,
        'dodec': 12,
        'tridec': 13,
        'tetradec': 14,
        'pentadec': 15,
        'hexadec': 16,
        'heptadec': 17,
        'octadec': 18,
        'nonadec': 19,
        'eicos': 20,
        'icos': 20,
        'heneicos': 21,
        'docos': 22,
        'tricos': 23,
        'tetracos': 24,
        'pentacos': 25,
        'hexacos': 26,
        'heptacos': 27,
        'octacos': 28,
        'nonacos': 29,
        'triacont': 30
    },
    'uns': {
        '': 1,
        'adi': 2,
        'atri': 3,
        'atetra': 4,
        'apenta': 5,
        'ahexa': 6,
        'ahepta': 7,
        'aocta': 8
    },
    'end': {
        'enoate': 1,
        'anoate': 0,
        'enoic acid': 1,
        'anoic acid': 0
    }
}

fa_greek = {}
for cc, uns, end in itertools.product(
    fa_greek_parts['cc'].items(),
    fa_greek_parts['uns'].items(),
    fa_greek_parts['end'].items()):
    
    if len(uns[0]) and end[1] == 0:
        continue
    
    fa_greek['%s%s%s' % (cc[0], uns[0], end[0])] = (cc[1], uns[1] * end[1])

