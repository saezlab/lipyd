#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `ltp` python module
#
#  Copyright (c) 2018 - EMBL
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
from future.utils import iteritems
from past.builtins import xrange, range, reduce

import os
import time
import re
import imp
import sys
import copy
import itertools
import collections

import ltp.table as table
import ltp.settings as ltpsettings
import lipyd.common as common


class ResultsReprocessor(object):
    """
    Reprocesses a directory with top features xlsx tables.
    Finds mgf files for each protein for all protein containing fractions.
    Opens the xlsx files, does custom operations, writes data in new columns
    and saves xlsx files in a new directory with preserving formatting.
    """
    
    rexlsx = re.compile(r'([A-Z0-9]+)_top_features_.*xlsx')
    remgf  = re.compile(r'([A-Z0-9]+)_(pos|neg)_([A-Z])([0-9]{1,2})\.mgf')
    refr   = re.compile(r'([A-Z])([0-9]{1,2})')
    
    def __init__(
            self,
            source_dir,
            mgfdir,
            target_dir = None,
            screen = 'invitro',
            reprocessor = 'MS2Reprocessor',
        ):
        
        self.screen = screen
        self.source_dir = source_dir
        self.mgfdir = mgfdir
        self.target_dir = target_dir or '%s__%s' % (
            target_dir,
            time.strftime('%Y.%m.%d_%H.%M')
        )
        
        if not os.path.exists(self.target_dir):
            
            os.mkdir(self.target_dir)
        
        self.fractionsf = (
            ltpsettings.get('protein_containing_fractions_%s' % self.screen)
        )
        
        self.reprocessor = (
            getattr(table, reprocessor)
                if type(reprocessor) is common.basestring else
            reprocessor
        )
    
    def reload(self):
        
        modname = self.__class__.__module__
        mod = __import__(modname, fromlist = [modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)
    
    def main(self):
        
        self.read_protein_containing_fractions()
        self.collect_files()
        self.collect_mgfs()
    
    def read_protein_containing_fractions(self):
        
        result = collections.defaultdict(dict)
        
        def fraction_tuple(fr_str):
            
            fr = self.refr.match(fr_str).groups()
            
            return fr[0], int(fr[1])
        
        with open(self.fractionsf, 'r') as fp:
            
            hdr = fp.readline().split('\t')
            
            for l in fp:
                
                l = l.split('\t')
                
                if self.screen == 'invitro':
                    
                    fr = fraction_tuple(l[1])
                    
                    result[l[0]][fr] = int(l[2])
                
                if self.screen == 'invivo':
                    
                    for i, val in enumerate(l[1:]):
                        
                        if val != 'None':
                            
                            fr = fraction_tuple(hdr[i + 1].upper())
                            
                            result[l[0]][fr] = int(val)
        
        self.fractions = result
    
    def collect_files(self):
        
        result = {}
        
        for fname in os.listdir(self.source_dir):
            
            protein = self.rexlsx.match(fname)
            
            if not protein:
                
                continue
            
            protein = protein.groups()[0]
            
            result[protein] = fname
        
        self.sources = result
    
    def collect_mgfs(self):
        
        result = {}
        
        for mgf in os.listdir(self.mgfdir):
            
            match = self.remgf.search(mgf)
            
            if not match:
                
                continue
            
            protein, ionmode, fracrow, fraccol = match.groups()
            
            result[(protein, ionmode, (fracrow, fraccol))] = mgf
        
        self.mgfs = result
    
    def iter_tables(self):
        
        for protein, xlsx in iteritems(self.xlsx):
            
            xlsxpath = os.path.join(self.source_dir, xlsx)
            
            for ionmode in ('neg', 'pos'):
                
                mgfpaths = dict(
                    (
                        (ionmode, fr),
                        os.path.join(
                            self.mgfdir,
                            self.mgfs[(protein, ionmode, fr)],
                        )
                    )
                    for fr, pcont in iteritems(self.fractions[protein])
                    if pcont
                )
                
                yield xlsxpath, mgfpaths
    
    def reprocess(self):
        
        for xlsxpath, mgfpaths in self.iter_tables():
            
            reproc = self.reprocessor(
                infile = xlsxpath,
                outdir = self.target_dir,
                mgf_files = mgfpaths,
            )
            
            reproc.main()


#TODO: make sure all code below works and remove it from lipyd.main.Screening
class Fractions(object):
    
    refr = re.compile(r'([A-Z])([0-9]{1,2})-?([A-Z]?)([0-9]{0,2})')
    
    def __init__(
            self,
            screen,
        ):
        
        self.screen = screen
        
        # settings
        self.fractionsf = ltpsettings.get('%s_fractionsf' % self.screen)
        self.pcont_fracs_from_abs = (
            ltpsettings.get('pcont_fracs_from_abs')[self.screen]
        )
        self.use_manual_ppratios = ltpsettings.get('use_manual_ppratios')
        self.manual_ppratios_xls = ltpsettings.get('manual_ppratios_xls')
        self.manual_ppratios_xls_cols = (
            ltpsettings.get('manual_ppratios_xls_cols')
        )
    
    def reload(self):
        
        modname = self.__class__.__module__
        mod = __import__(modname, fromlist = [modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)
    
    def read_fractions(self):
        
        self.read_protein_containing_fractions()
        self.upper_fractions()
        self.raw_sec_absorbances()
        self.absorbances_by_fractions()
        if self.pp_do_correction:
            self.pp_baseline_correction()
            self.pp_background()
            self.pp_background_correction()
        else:
            self.abs_by_frac_c = self.abs_by_frac
        self.mean_pp()
        if self.pcont_fracs_from_abs:
            self.protein_containing_fractions_from_absorbances()
        if not hasattr(self, 'fractions_upper'):
            self.fractions_upper = self.fractions
    
    def read_protein_containing_fractions(self):
        """
        Reads from file the sample/control annotations for each proteins.
        Also reads the highest fractions which is used if we use
        the profile evaluation method from Enric.
        """
        
        def get_names(names):
            
            return \
                filter(
                    lambda n:
                        n.upper() == n,
                    map(
                        lambda n:
                            n.replace(')', '').strip(),
                        names.split('(')
                    )
                )
        
        def get_fractions(frlist):
            
            return \
                list(
                    itertools.chain(
                        *map(
                            lambda i:
                                map(
                                    lambda num:
                                        '%s%s' % (i[0], num),
                                    xrange(int(i[1]), int(i[3]) + 1)
                                ) if len(i[2]) and len(i[3]) else \
                                ['%s%s' % (i[0], i[1])],
                            self.refr.findall(frlist)
                        )
                    )
                )
        
        def get_hpeak(hfracs):
            
            return (
                list(
                    map(
                        lambda hf:
                            hf.split('-'),
                        hfracs.split(':')
                    )
                )
            )
        
        # ##
        
        if not self.pcont_fracs_from_abs:
            
            data = {}
            
            with open(self.fractionsf, 'r') as fp:
                
                hdr = fp.readline().strip()
                
                col2fr = dict(
                    map(
                        lambda fr:
                            (fr[0], fr[1].upper()),
                        enumerate(hdr.split(',')[1:])
                    )
                )
                
                for l in fp:
                    
                    l = l.strip().split(',')
                    
                    data[l[0].replace('"', '').upper()] = (
                        dict(
                            map(
                                lambda x:
                                    (col2fr[x[0]], common.to_int(x[1])),
                                filter(
                                    lambda x:
                                        x[1] != '',
                                    enumerate(l[1:])
                                )
                            )
                        )
                    )
            
            self.fractions = data
            
        else:
            
            if self.fractionsf.endswith('xlsx'):
                
                tab = common.read_xls(self.fractionsf, sheet = 'status')
                
                self._fractions, self.hpeak = (
                    tuple(
                        map(
                            dict,
                            zip(
                                *itertools.chain(
                                    *map(
                                        lambda l:
                                            list(
                                                map(
                                                    lambda n:
                                                        [
                                                            (
                                                                n, # protein name
                                                                get_fractions(l[15])
                                                            ),
                                                            (
                                                                n, # protein name
                                                                get_hpeak(l[18])
                                                            )
                                                        ],
                                                    get_names(l[0])
                                                )
                                            ),
                                        filter(
                                            lambda l:
                                                len(l[15]) and \
                                                    l[15].strip() != 'NA' and (
                                                        l[15][0].upper() ==
                                                        l[15][0]) and (
                                                    l[15][1].isdigit()),
                                            tab[1:]
                                        )
                                    )
                                )
                            )
                        )
                    )
                )
                
            else:
                
                with open(self.fractionsf, 'r') as fp:
                    
                    self._fractions = \
                        dict(
                            map(
                                lambda l:
                                    (l[0], l[1:]),
                                map(
                                    lambda l:
                                        l.split(';'),
                                    filter(
                                        len,
                                        fp.read().split('\n')
                                    )
                                )
                            )
                        )
    
    def upper_fractions(self):
        """
        Creates the dict ``fractions_upper`` of fractions with
        uppercase LTP names, with same content as ``fractions``.
        """
        
        if hasattr(self, 'fractions'):
            
            self.fractions_upper = dict(
                (l.upper(), s)
                for l, s in
                iteritems(self.fractions)
            )
    
    def protein_containing_fractions_from_absorbances(self, abscol = 1):
        """
        Sets the protein containing fractions (`fractions`) based on
        the UV absorbances, considering protein containing those whit a
        value higher than the maximum / minratio.
        """
        
        self.fractions = dict(
            map(
                lambda protein: (protein, {}),
                self.abs_by_frac_c.keys()
            )
        )
        
        for protein, ab in iteritems(self.abs_by_frac_c):
            pr = sorted(list(ab.values())[0].keys()) # the UV profile IDs,
                                                     # e.g. 3 = 260 nm
            abs_max = \
                max(
                    map(
                        lambda a:
                            np.mean(a[abscol]),
                        ab.values()
                    )
                )
            
            self.fractions[protein] = \
                dict(
                    map(
                        lambda f:
                            (
                                f[0],
                                int(np.mean(f[1][abscol]) > \
                                    abs_max / self.pp_minratio)
                            ),
                        iteritems(ab)
                    )
                )
        
        for protein, frs in iteritems(self.fractions):
            if protein in self._fractions:
                for fr, val in iteritems(frs):
                    if fr not in self._fractions[protein]:
                        self.fractions[protein][fr] = 0
                    else:
                        self.fractions[protein][fr] = 1
    
    def fractions_marco(self):
        """
        Reads manually defined protein peak ratios and protein containing
        fractions.
        """
        
        self.fractions_original = \
            copy.deepcopy(self.fractions_upper
                          if hasattr(self, 'fractions_upper')
                          else self.fractions)
        
        if not self.use_manual_ppratios:
            return None
        
        refrac = re.compile(r'([AB])([0-9]{1,2})')
        tbl = self.read_xls(self.manual_ppratios_xls)[1:]
        ix = self.manual_ppratios_xls_cols
        for l in tbl:
            protein = l[ix[0]].split('=')[0].strip()
            if not len(protein):
                continue
            fracs = set(map(
                lambda fr:
                    '%s%u' % (fr[0].lower(), int(fr[1])),
                refrac.findall(l[ix[1]])
            ))
            for i, fr in enumerate(self.fracs):
                if fr in fracs:
                    if self.fractions_upper[protein][i + 1] != 1:
                        sys.stdout.write('\t:: Setting fraction %s at %s to 1, '\
                            'this was %s before\n' % \
                            (fr, protein, self.fractions_upper[protein][i + 1]))
                        self.fractions_upper[protein][i + 1] = 1
                elif self.fractions_upper[protein][i + 1] == 1:
                    sys.stdout.write('\t:: Setting fraction %s at %s to 0, '\
                        'this was %s before\n' % \
                        (fr, protein, self.fractions_upper[protein][i + 1]))
                    self.fractions_upper[protein][i + 1] = 0
