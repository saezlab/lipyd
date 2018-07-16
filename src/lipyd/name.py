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

from __future__ import print_function
from past.builtins import xrange, range, reduce

from future.utils import iteritems

import sys
import imp
import re
import itertools

import lipyd.settings as settings


class LipidNameProcessor(object):
    
    def __init__(
            self,
            database = 'swisslipids',
            with_alcohols = True,
            with_coa = True,
            iso = False
        ):
        """
        Processes lipid names used in databases. Converts names to the
        standard used in this module and extracts carbon count and
        unsaturation information and other features.
        """
        
        self.database = database.lower()
        self.with_alcohols = with_alcohols
        self.with_coa = with_coa
        self.iso = iso
        self.lipnamesf = settings.get('lipnamesf')
        self.adducts_constraints = settings.get('adducts_constraints')
        self.recount1 = re.compile(
            r'\(([POdt]?)-?([0-9]{1,2}):([0-9]{1,2})\)'
        )
        self.recount2 = re.compile(
            r'\(([POdt]?)-?([0-9]{1,2}):([0-9]{1,2})/'
            r'([POdt]?)-?([0-9]{1,2}):([0-9]{1,2})/?'
            r'([POdt]?)-?([0-9]{0,2}):?([0-9]{0,2})\)'
        )
        self.recount3 = re.compile(
            r'\(([POdt]?)-?([0-9]{1,2}):([0-9]{1,2})[/_]?'
            r'([POdt]?)-?([0-9]{0,2}):?([0-9]{0,2})[/_]?'
            r'([POdt]?)-?([0-9]{0,2}):?([0-9]{0,2})[/_]?'
            r'([POdt]?)-?([0-9]{0,2}):?([0-9]{0,2})\)'
        )
        self.recount4 = re.compile(
            r'\(?'
            r'((?:[0-9]+-)?[POdt]?)-?([0-9]{1,2}):([0-9]{1,2})'
            r'\(?([0-9EZ,]*)\)?((?:[-\(]2OH\)?)?)[/_]?'
            r'((?:[0-9]+-)?[POdt]?)-?([0-9]{0,2}):?([0-9]{0,2})'
            r'\(?([0-9EZ,]*)\)?((?:[-\(]2OH\)?)?)[/_]?'
            r'((?:[0-9]+-)?[POdt]?)-?([0-9]{0,2}):?([0-9]{0,2})'
            r'\(?([0-9EZ,]*)\)?((?:[-\(]2OH\)?)?)[/_]?'
            r'((?:[0-9]+-)?[POdt]?)-?([0-9]{0,2}):?([0-9]{0,2})'
            r'\(?([0-9EZ,]*)\)?((?:[-\(]2OH\)?)?)\)?'
        )
        self.reme     = re.compile(r'methyl|ethyl')
        self.rebr2    = re.compile(
            r'(1(?:,2-di)?)-\(((?:[0-9]{0,2}[-]?methyl|ethyl)?)[A-z0-9-]+\)'
            r'-([2,3]{1,3}(?:-di)?)-'
            r'\(((?:[0-9]{0,2}[-]?methyl|ethyl)?)[A-z0-9-]+\)'
        )
        self.gen_fa_greek()
        self.read_lipid_names()
    
    def reload(self, children = False):
        modname = self.__class__.__module__
        mod = __import__(modname, fromlist=[modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)
    
    def read_lipid_names(self, add_fa = True):
        """
        Reads annotations for lipid classes:
        - full names
        - short notations
        - database keywords
        (to process long names from SwissLipids and LipidMaps)
        - most abundant adducts
        
        The input file is given by the ``lipnamesf`` attribute.
        """
        result = {}
        
        with open(self.lipnamesf, 'r') as f:
            nul = f.readline()
            for l in f:
                l = l.strip().split('\t')
                result[l[0]] = {
                    'full_name': l[1],
                    'swl': self.process_db_keywords(l[2]),
                    'lmp': self.process_db_keywords(l[3]),
                    'pos_adduct': (
                            l[4]
                        if l[4] != 'ND' and self.adducts_constraints else
                            None
                        ),
                    'neg_adduct': (
                            l[5]
                        if l[5] != 'ND' and self.adducts_constraints else
                            None
                        )
                }
        
        result['FA'] = {'full_name': 'Fatty acid', 'swl': [], 'lmp': [],
                        'pos_adduct': None, 'neg_adduct': None}
        
        self.lipnames = result
    
    @staticmethod
    def process_db_keywords(kwdstr):
        
        return \
            list(
                map(
                    lambda kwdset:
                        {
                            'neg':
                                list(
                                    map(
                                        lambda kwd:
                                            kwd[1:],
                                        list(
                                            filter(
                                                lambda kwd:
                                                    len(kwd) and \
                                                    kwd[0] == '!',
                                                kwdset.split(';')
                                            )
                                        )
                                    )
                                ),
                            'pos':
                                list(
                                    filter(
                                        lambda kwd:
                                            len(kwd) and kwd[0] != '!',
                                        kwdset.split(';')
                                    )
                                )
                        },
                    kwdstr.split('|')
                )
            )
    
    @staticmethod
    def prefix_proc(prfx):
        
        return 'O-' if prfx in {'O', 'P'} else prfx
    
    def carbon_counts(self, name, ccexp = 2):
        """
        Processes carbon and unsaturation counts from name.
        
        Args
        ----
        :param str name:
            Lipid name.
        :param int ccexp:
            Expected number of fatty acyl or other residues constaining
            aliphatic chain. E.g. for DAG this should be 2 and for TAG 3
            as the latter has 3 fatty acyls.
        """
        
        # regex finds the total carbon count
        cc1 = self.recount1.findall(name)
        # regex finds 2-3 fatty acids
        cc3 = self.recount3.findall(name)
        
        # the total carbon count
        ccpart = (
            [self.prefix_proc(cc1[0][0]), int(cc1[0][1]), int(cc1[0][2])]
            if len(cc1) else
            None
        )
        
        faccparts = []
        
        if ccexp and cc3 and cc3[0][(ccexp - 1) * 3 + 1]:
            
            for i in range(ccexp):
                
                if cc3[0][i * 3 + 1]:
                    
                    faccparts.append((
                            self.prefix_proc(cc3[0][i * 3]),
                            int(cc3[0][i * 3 + 1]),
                            int(cc3[0][i * 3 + 2])
                    ))
            
            if len(faccparts) != ccexp:
                
                faccparts = None
        
        return ccpart, faccparts
    
    def isomeric_carbon_counts(self, name):
        
        icc = self.recount4.findall(name)
        
        if icc:
            
            try:
                
                icc = [
                    (
                        # the usual prefix, carbon, unsat tuple:
                        (
                            self.prefix_proc(icc[0][i]),
                            int(icc[0][i + 1]),
                            int(icc[0][i + 2])
                        ),
                        # and the isomeric conformation string:
                        icc[0][i + 3],
                        # and the hydroxyl group on C2:
                        bool(icc[0][i + 4])
                    )
                    for i in xrange(0, 16, 5)
                    # keep only existing aliphatic chains
                    if icc[0][i + 1] and icc[0][i + 2]
                ]
            
            except Exception as e:
                
                sys.stdout.write(''.join([
                    '\n\n\n',
                    '!!! Incomprehensible results at isomeric carbon counts:',
                    '\nRegex match: %s\nName processed: %s\n' % (
                        icc, name
                    ),
                    '\n\n\n'
                ]))
        
        return icc
    
    def headgroup_from_lipid_name(self, name, database = None):
        """
        For one database record attempts to identify the lipid class
        by looking up keywords.
        Calls greek name identification, greek fatty acid names are
        identified as 'FA'.
        """
        
        database = database or self.database
        
        db = 'lmp' if database.lower() == 'lipidmaps' else 'swl'
        for shortname, spec in iteritems(self.lipnames):
            for kwset in spec[db]:
                matched = [kw in name for kw in kwset['pos']]
                if sum(matched) == len(kwset['pos']) and \
                    sum(matched) > 0:
                    matched = [kw in name for kw in kwset['neg']]
                    if sum(matched) == 0:
                        return shortname
        
        return self.process_fa_name(name)
    
    def process_fa_name(self, name):
        """
        Identifies fatty acids based on their greek name.
        """
        
        return (
            'FA'
                if name in self.fa_greek
                or 'FA' in name
                or 'atty acid' in name
            else 'FAL'
                if self.with_alcohols and (
                    name in self.fal_greek or
                    'atty alcohol' in name
                )
            else 'FACoA'
                # TODO: capture carbon counts of acyl-CoAs
                if self.with_coa and (
                    name in self.facoa_greek or
                    'oyl-CoA' in name
                )
            else None
        )
    
    def fa_greek_cc(self, name):
        
        cc1, cc2, cc3 = None, [], []
        
        try:
            
            name1 = name.split('-')[1] if '-' in name else name
            
            for dct in ['fa', 'fal', 'facoa']:
                
                if name1 in getattr(self, '%s_greek' % dct):
                    
                    cc1 = getattr(self, '%s_greek' % dct)[name1]
                    cc1 = ('', cc1[0], cc1[1])
                    cc2 = [cc1]
                    cc3 = [(
                        cc1,
                        name.split(')')[0][1:] if '(' in name else '',
                        False
                    )]
        
        except:
            pass
        
        return cc1, cc2, cc3
    
    def test_branched(self, name):
        """
        Tells if a lipid might contain branched aliphatic chains.
        """
        
        return bool(self.reme.search(name))
    
    def process(self, name, database = None):
        """
        The main method of this class. Processes a lipid name string
        and returns a standard name, prefix, carbon counts and
        unsaturations.
        
        Args
        ----
        :param str name:
            One or more names to process. Single result will be returned
            and names will be attempted to be processed one after the other
            until processing is successful. Names in one string can be
            separated by `|`.
        """
        
        database = database or self.database
        
        hg, cc1, cc2, icc = None, None, None, None
        
        hg = self.headgroup_from_lipid_name(name, database = database)
        
        # try greek fatty acyl carbon counts:
        
        
        if not hg and self.iso and database == 'swisslipids':
            
            try:
                
                for name0 in name.split('|'):
                    
                    fa_greek = name0.split('-')[1]
                    hg = self.process_fa_name(fa_greek)
                    if hg:
                        
                        break
                
            except:
                
                pass
        
        for n in name.split('|'):
            
            lyso =  hg and hg[:4] == 'Lyso'
            
            # how many aliphatic chains this molecule has
            ccexp = (
                    1
                if hg in {'FA', 'MAG'} or lyso else
                    3
                if hg == 'TAG' else
                    2
            )
            
            if lyso and '0:0' in n:
                # SwissLipids adds `0:0` to lyso glycerolipids
                ccexp += 1
            
            if hg == 'BMP' and '0:0' in n:
                # SwissLipids shows 4 acyl residues for BMP
                # 2 of them are `0:0`
                ccexp = 4
            
            _cc1, _cc2 = self.carbon_counts(n, ccexp = ccexp)
            
            if self.iso and not icc:
                
                icc = self.isomeric_carbon_counts(n)
                
                if icc and not cc2:
                    
                    cc2 = [i[0] for i in icc]
            
            cc1 = cc1 or _cc1
            cc2 = cc2 or _cc2
            
            if cc2 and cc1 and (not self.iso or icc):
                
                break
        
        if not cc1 and cc2:
            
            cc1 = [
                ''.join(i[0] for i in cc2),
                sum(i[1] for i in cc2),
                sum(i[2] for i in cc2)
            ]
        
        if hg in set(['FA', 'FAL', 'FACoA']) and not cc1:
            
            for name0 in name.split('|'):
                
                cc1, cc2, icc = self.fa_greek_cc(name0)
                
                if cc1:
                    break
        
        if cc1 and cc1[0] == 'd' and cc1[2] == 0:
            
            cc1[0] = 'DH'
        
        # if the sphingolipid headgroup does not contain the prefix:
        if hg and cc1 and cc1[0]:
            
            for prfx in ('d', 't', 'DH'):
                
                if cc1[0] == prfx and not hg.startswith(prfx):
                    
                    hg = '%s%s' % (prfx, hg)
                    break
        
        return hg, cc1, cc2, icc
    
    def gen_fa_greek(self):
        """
        Generates a list of greek fatty acid names with their
        carbon counts and unsaturations.
        """
        
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
        
        fal_greek_end = {}
        fal_greek_end['anol'] = 0
        fal_greek_end['enol'] = 1
        
        facoa_greek_end = {}
        facoa_greek_end['anoyl-CoA'] = 0
        facoa_greek_end['enoyl-CoA'] = 1
        
        self.fa_greek  = {}
        self.fal_greek = {}
        self.facoa_greek = {}
        
        for cc, uns, end in itertools.product(
            fa_greek_parts['cc'].items(),
            fa_greek_parts['uns'].items(),
            fa_greek_parts['end'].items()):
            
            if len(uns[0]) and end[1] == 0:
                continue
            
            self.fa_greek['%s%s%s' % (cc[0], uns[0], end[0])] = (cc[1], uns[1] * end[1])
        
        for cc, uns, end in itertools.product(
            fa_greek_parts['cc'].items(),
            fa_greek_parts['uns'].items(),
            fal_greek_end.items()):
            
            if len(uns[0]) and end[1] == 0:
                continue
            
            self.fal_greek['%s%s%s' % (cc[0], uns[0], end[0])] = (cc[1], uns[1] * end[1])
        
        for cc, uns, end in itertools.product(
            fa_greek_parts['cc'].items(),
            fa_greek_parts['uns'].items(),
            facoa_greek_end.items()):
            
            if len(uns[0]) and end[1] == 0:
                continue
            
            self.facoa_greek['%s%s%s' % (cc[0], uns[0], end[0])] = (cc[1], uns[1] * end[1])
