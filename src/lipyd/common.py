#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `lipyd` python module
#
#  Copyright (c) 2015-2017 - EMBL
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
import os
import warnings
import traceback
import itertools
import xlrd
import openpyxl

ROOT = os.path.abspath(os.path.dirname(__file__))

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


class _const:
    
    class ConstError(TypeError):
        
        pass

    def __setattr__(self, name, value):
        
        if name in self.__dict__:
            
            raise(self.ConstError, "Can't rebind const(%s)" % name)
        
        self.__dict__[name] = value


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

count_prefix = {
    1: 'mono',
    2: 'di',
    3: 'tri',
    4: 'tetra',
    5: 'penta',
    6: 'hexa',
    7: 'hepta',
    8: 'octa',
    9: 'nona'
}

fa_greek = {}
for cc, uns, end in itertools.product(
    fa_greek_parts['cc'].items(),
    fa_greek_parts['uns'].items(),
    fa_greek_parts['end'].items()):
    
    if len(uns[0]) and end[1] == 0:
        continue
    
    fa_greek['%s%s%s' % (cc[0], uns[0], end[0])] = (cc[1], uns[1] * end[1])


adducts = {
    'pos': ['[M+H]+', '[M+NH4]+', '[M+Na]+'],
    'neg': ['[M-H]-', '[M+HCOO]-']
}


ad2ex = {
    1: {
        'pos': {
            '[M+H]+': 'remove_h',
            '[M+NH4]+': 'remove_nh4',
            '[M+Na]+': 'remove_na',
        },
        'neg': {
            '[M-H]-': 'add_h',
            '[M+HCOO]-': 'remove_fo'
        }
    },
    2: {
        'pos': {},
        'neg': {
            '[M-2H]2-': 'add_2h'
        }
    },
    3: {
        'pos': {},
        'neg': {
            '[M-3H]3-': 'add_3h'
        }
    }
}

exact_method = {
    '[M+OAc]-': 'remove_ac',
    '[M+H]+': 'remove_h',
    '[M-H]-': 'add_h',
    '[M+HCOO]-': 'remove_fo',
    '[M+NH4]+': 'remove_nh4',
    '[M+Na]+': 'remove_na'
}

# method names to convert between exact and adduct masses
adduct_method = {
    '[M+OAc]-': 'add_ac',
    '[M+H]+': 'add_h',
    '[M-H]-': 'remove_h',
    '[M+HCOO]-': 'add_fo',
    '[M+NH4]+': 'add_nh4',
    '[M+Na]+': 'add_na'
}

ex2ad = {
    1: {
        'pos': {
            '[M+H]+': 'add_h',
            '[M+NH4]+': 'add_nh4',
            '[M+Na]+': 'add_na'
        },
        'neg': {
            '[M-H]-': 'remove_h',
            '[M+HCOO]-': 'add_fo'
        }
    },
    2: {
        'pos': {},
        'neg': {
            '[M-2H]2-': 'remove_2h'
        }
    },
    3: {
        'pos': {},
        'neg': {
            '[M-3H]3-': 'remove_3h'
        }
    }
}


def read_xls(xls_file, sheet = 0, csv_file = None,
    return_table = True):
    """
    Generic function to read MS Excel XLS file, and convert one sheet
    to CSV, or return as a list of lists
    """
    table = []
    try:
        book = xlrd.open_workbook(xls_file, on_demand = True)
        try:
            if type(sheet) is int:
                sheet = book.sheet_by_index(sheet)
            else:
                sheet = book.sheet_by_name(sheet)
        except xlrd.biffh.XLRDError:
            sheet = book.sheet_by_index(0)
        table = [[unicode(c.value) \
            for c in sheet.row(i)] \
            for i in xrange(sheet.nrows)]
    except IOError:
        sys.stdout.write('No such file: %s\n' % xls_file)
        sys.stdout.flush()
    except:
        try:
            book = openpyxl.load_workbook(filename = xls_file,
                read_only = True)
        except:
            sys.stdout.write('\tCould not open xls: %s\n' % xls_file)
            if not os.path.exists(xls_file):
                sys.stdout.write('\tFile does not exist.\n')
            sys.stdout.flush()
        try:
            if type(sheet) is int:
                sheet = book.worksheets[sheet]
            else:
                sheet = book[sheet]
        except:
            sheet = book.worksheets[0]
        cells = sheet.get_squared_range(1, 1,
            sheet.max_column, sheet.max_row)
        table = map(lambda row:
            map(lambda c:
                unicode(c.value) if c.value else '',
                row
            ),
            cells
        )
    if csv_file:
        with open(csv_file, 'w') as csv:
            csv.write('\n'.join(['\t'.join(r) for r in table]))
    if not return_table:
        table = None
    if 'book' in locals() and hasattr(book, 'release_resources'):
        book.release_resources()
    return table
