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

import os
import sys
import imp
import re
import copy
import struct
import itertools
import collections

import numpy as np
import pandas as pd

try:
    import pybel
except:
    sys.stdout.write('\t:: No module `pybel` available.\n')

import lipyd._curl as _curl
import lipyd.common as common
import lipyd.settings as settings
import lipyd.mz as mzmod
import lipyd.progress as progress
import lipyd.sdf as sdf
import lipyd.lipid as lipid
import lipyd.lookup as lookup


class Reader(object):
    
    def __init__(self):
        
        self.load()
        self.process()
        return self.__iter__()
    
    def reload(self, children = False):
        
        modname = self.__class__.__module__
        mod = __import__(modname, fromlist=[modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)
    
    def iterrows(self):
        
        for mol in self:
            
            yield [
                mol[2].data['LM_ID'],
                'Species',
                '|'.join(self.names(mol[2])),
                mol[2].data['INCHI_KEY'],
                '',
                mol[0],
                mol[2].data['PUBCHEM_CID']
            ]


class LipidMapsOld(object):
    
    def __init__(self):
        
        self.url = common.get_param('lipidmaps_url')
        self.lipidmaps_fname  = common.get_param('lipidmaps_fname')
        
        SdfReader.__init__(self, name = 'LipidMaps')
        Reader.__init__(self)
    
    @staticmethod
    def names(mol):
        
        return (mol[k]
                for k in
                ('COMMON_NAME', 'SYNONYMS', 'SYSTEMATIC_NAME')
                if k in mol)
    
    def load(self):
        
        c = _curl.Curl(self.url, large = True,
                       silent = False, files_needed = [self.lipidmaps_fname])
        fpa = c.result[fn]
        self.fname = os.path.join('cache', fn)
        
        with open(self.fname, 'wb') as fpe:
            
            while True:
                
                block = fpa.read(4096)
                if not block:
                    break
                
                fpe.write(block)
        
        c.close()
    
    def iterrows(self):
        
        for mol in self:
            
            yield [
                    mol[2].data['LM_ID'],
                    'Species',
                    '|'.join(self.names(mol[2])),
                    mol[2].data['INCHI_KEY'],
                    '',
                    mol[0],
                    mol[2].data['PUBCHEM_CID']
                ]


class LipidMaps(sdf.SdfReader):
    
    def __init__(self, extract_file = True):
        
        self.url   = settings.get('lipidmaps_url')
        self.fname = settings.get('lipidmaps_fname')
        self.curl  = _curl.Curl(self.url, large = True, silent = False)
        
        if extract_file:
            
            self.efname = os.path.join('cache', self.fname.split('/')[-1])
            
            with open(self.efname, 'wb') as efp:
                
                for l in self.curl.result[self.fname]:
                    
                    efp.write(l)
            
            efp = open(os.path.join('cache', self.fname.split('/')[-1]), 'rb')
            sdf.SdfReader.__init__(self, efp)
        
        else:
            sdf.SdfReader.__init__(self, self.curl.result[self.fname])


class SwissLipids(Reader):
    
    def __init__(self, levels = set(['Species']), silent = False,
                 nameproc_args = {}, branched = False):
        
        self.silent = silent
        self.nameproc_args = nameproc_args
        self.set_levels(levels)
        self.url = settings.get('swisslipids_url')
        self.load()
        self.make_index()
    
    def set_levels(self, levels):
        """
        Sets levels considered. Levels in SwissLipids are `Species`,
        `Molecular subspecies`, `Structural subspecies` and
        `Isomeric subspecies`.
        
        :param set levels: A set of one or more of the levels above.
        """
        
        self.levels = levels
        self.init_name_processor()
    
    def init_name_processor(self):
        """
        Creates a `LipidNameProcessor` instance to process lipid names
        at indexing.
        """
        
        self.nameproc = LipidNameProcessor(
            iso = 'Isomeric subspecies' in self.levels,
            **self.nameproc_args
        )
    
    def load(self):
        
        self.close_gzfile()
        
        self._curl = _curl.Curl(self.url, silent = False,
                                compr = 'gz', large = True)
        self._gzfile = self._curl.result
    
    def iterfile(self):
        
        self._plainfile.seek(0)
        _ = self._plainfile.readline()
        
        for line in self._plainfile:
            
            yield line
    
    @staticmethod
    def names(line):
        
        return '|'.join(line[2:5])
    
    def make_index(self):
        
        def cc2str(cc):
            
            return (
                '%s%s%u:%u' % (
                    cc[0],
                    '-' if cc[0] in {'O', 'P'} else '',
                    cc[1],
                    cc[2]
                )
            )
        
        self.close_plainfile()
        
        self.load()
        self.index = collections.defaultdict(lambda: set([]))
        self.hg_index      = collections.defaultdict(lambda: set([]))
        self.species_index = collections.defaultdict(lambda: set([]))
        self.subspec_index = collections.defaultdict(lambda: set([]))
        self.isomer_index  = collections.defaultdict(lambda: set([]))
        
        if not self.silent:
            
            self._gzfile.fileobj.seek(-4, 2)
            ucsize = struct.unpack('I', self._gzfile.fileobj.read(4))[0]
            self.prg = progress.Progress(ucsize, 'Indexing SwissLipids', 101)
        
        self._gzfile.fileobj.seek(0)
        
        self._plainfilename = '%s.extracted' % self._gzfile.name
        
        with open(self._plainfilename, 'wb') as fpp:
            
            offset = self._gzfile.tell()
            
            for l in self._gzfile:
                
                if not self.silent:
                    self.prg.step(len(l))
                
                ll = l.decode('ascii').split('\t')
                
                if ll[1] in self.levels:
                    
                    names = self.names(ll)
                    self.index[ll[0]].add(offset) # SwissLipids ID
                    self.index[ll[8]].add(offset) # SMILES
                    self.index[ll[10]].add(offset) # InChI key
                    
                    for n in names.split('|'):
                        
                        self.index[n].add(offset)
                    
                    hg, cc1, cc2, icc = self.nameproc.process(names)
                    
                    if hg:
                        
                        self.hg_index[hg].add(offset)
                    
                    if cc1:
                        
                        self.species_index[
                            '%s(%s)' % (hg, cc2str(cc1))
                        ].add(offset)
                    
                    if cc2:
                        
                        self.subspec_index[
                            '%s(%s)' % (
                                hg,
                                '/'.join(cc2str(a) for a in cc2)
                            )
                        ].add(offset)
                    
                    if icc:
                        
                        self.isomer_index[
                            '%s(%s)' % (
                                hg,
                                '/'.join('%s(%s)%s' % (
                                    cc2str(a[0]),
                                    a[1],
                                    '-2OH' if a[2] else ''
                                ) for a in icc)
                            )
                        ].add(offset)
                
                offset = self._gzfile.tell()
                fpp.write(l)
        
        if not self.silent:
            self.prg.terminate()
        
        self.index = dict(self.index)
        self._plainfile = open(self._plainfilename, 'r')
    
    def get_hg(self, hg):
        
        return self.get_record(hg, index = 'hg')
    
    def get_hg_obmol(self, hg):
        
        return self.get_obmol(hg, index = 'hg')
    
    def get_species(self, name):
        
        return self.get_record(name, index = 'species')
    
    def get_subspec(self, name):
        
        return self.get_record(name, index = 'subspec')
    
    def get_isomer(self, name):
        
        return self.get_record(name, index = 'isomer')
    
    def get_hg_obmol(self, hg):
        
        return self.get_obmol(hg, index = 'hg')
    
    def get_species_obmol(self, name):
        
        return self.get_obmol(name, index = 'species')
    
    def get_subspec_obmol(self, name):
        
        return self.get_obmol(name, index = 'subspec')
    
    def get_isomer_obmol(self, name):
        
        return self.get_obmol(name, index = 'isomer')
    
    def get_record(self, name, index = ''):
        
        indexname = '%s%sindex' % (index, '_' if index else '')
        index = getattr(self, indexname)
        
        if name in index:
            
            for offset in index[name]:
                
                self._plainfile.seek(offset)
                yield self._plainfile.readline().strip().split('\t')
    
    def get_obmol(self, name, index = ''):
        
        return [
            self.to_obmol(rec)
            for rec in self.get_record(name, index = index)
        ]
    
    @staticmethod
    def to_obmol(record):
        
        if not record[8]:
            
            return None
        
        mol = pybel.readstring('inchi', record[9])
        mol.db_id = record[0]
        mol.name  = record[3]
        mol.title = '|'.join(record[2:5])
        mol.chebi = record[24] if len(record) > 24 else ''
        mol.lipidmaps = record[25] if len(record) > 25 else ''
        mol.hmdb = record[26] if len(record) > 26 else ''
        mol.smiles = record[8]
        mol.swl_exact_mass = float(record[14]) if record[14] else None
        mol.swl_formula = record[11]
        mol.inchi = record[9]
        mol.inchikey = record[10]
        mol.level = record[1]
        
        return mol
    
    def __iter__(self):
        
        nosmiles = 0
        
        for line in self.iterfile():
            
            line = line.strip().split('\t')
            
            if len(line) > 22 and line[1] in self.levels:
                
                if not line[8]:
                    nosmiles += 1
                    continue
                
                mol = self.to_obmol(line)
                
                yield mol
    
    def __del__(self):
        
        self.close_gzfile()
        self.close_plainfile()
    
    def close_gzfile(self):
        
        if hasattr(self, '_gzfile'):
            
            self._gzfile.close()
    
    def close_plainfile(self):
        
        if hasattr(self, '_plainfile'):
            
            self._plainfile.close()
    
    def export_names(self, proc):
        
        with open('names.tmp', 'w') as fp:
            
            for i in self.__iter__():
                
                n = proc.process(i.title)
                
                fp.write('%s\t%s\n' % (i.title, str(n)))


class LipidNameProcessor(object):
    
    def __init__(self, with_alcohols = True, with_coa = True, iso = False):
        
        self.with_alcohols = with_alcohols
        self.with_coa = with_coa
        self.iso = iso
        self.lipnamesf = settings.get('lipnamesf')
        self.adducts_constraints = settings.get('adducts_constraints')
        self.recount1 = re.compile(r'\(([POdt]?)-?([0-9]{1,2}):([0-9]{1,2})\)')
        self.recount2 = re.compile(r'\(([POdt]?)-?([0-9]{1,2}):([0-9]{1,2})/'
                                   r'([POdt]?)-?([0-9]{1,2}):([0-9]{1,2})/?'
                                   r'([POdt]?)-?([0-9]{0,2}):?([0-9]{0,2})\)')
        self.recount3 = re.compile(r'\(([POdt]?)-?([0-9]{1,2}):([0-9]{1,2})[/_]?'
                                   r'([POdt]?)-?([0-9]{0,2}):?([0-9]{0,2})[/_]?'
                                   r'([POdt]?)-?([0-9]{0,2}):?([0-9]{0,2})[/_]?'
                                   r'([POdt]?)-?([0-9]{0,2}):?([0-9]{0,2})\)')
        self.recount4 = re.compile(r'\(?([POdt]?)-?([0-9]{1,2}):([0-9]{1,2})\(?([0-9EZ,]*)\)?((?:-2OH)?)[/_]?'
                                   r'([POdt]?)-?([0-9]{0,2}):?([0-9]{0,2})\(?([0-9EZ,]*)\)?((?:-2OH)?)[/_]?'
                                   r'([POdt]?)-?([0-9]{0,2}):?([0-9]{0,2})\(?([0-9EZ,]*)\)?((?:-2OH)?)[/_]?'
                                   r'([POdt]?)-?([0-9]{0,2}):?([0-9]{0,2})\(?([0-9EZ,]*)\)?((?:-2OH)?)\)?')
        self.reme     = re.compile(r'methyl|ethyl')
        self.rebr2    = re.compile(r'(1(?:,2-di)?)-\(((?:[0-9]{0,2}[-]?methyl|ethyl)?)[A-z0-9-]+\)'
                                   r'-([2,3]{1,3}(?:-di)?)-\(((?:[0-9]{0,2}[-]?methyl|ethyl)?)[A-z0-9-]+\)')
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
                    'pos_adduct': l[4] if l[4] != 'ND' and self.adducts_constraints else None,
                    'neg_adduct': l[5] if l[5] != 'ND' and self.adducts_constraints else None
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
    
    def carbon_counts(self, name, ccexp = 2):
        
        # regex finds the total carbon count
        cc1 = self.recount1.findall(name)
        # regex finds 2-3 fatty acids
        cc3 = self.recount3.findall(name)
        
        # the total carbon count
        ccpart = (
            [cc1[0][0], int(cc1[0][1]), int(cc1[0][2])]
            if len(cc1) else
            None
        )
        
        faccparts = []
        
        if ccexp and cc3 and cc3[0][(ccexp - 1) * 3 + 1]:
            
            for i in range(ccexp):
                
                if cc3[0][i * 3 + 1]:
                    
                    faccparts.append((
                            cc3[0][i * 3],
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
                        (icc[0][i], int(icc[0][i + 1]), int(icc[0][i + 2])),
                        # and the isomeric conformation string:
                        icc[0][i + 3],
                        # and the hydroxyl group on C2:
                        bool(icc[0][i + 4])
                    )
                    for i in xrange(0, 16, 5)
                    # keep only existing aliphatic chains
                    if icc[0][i + 1]
                ]
            
            except:
                
                sys.stdout.write(''.join([
                    '\n\n\n',
                    '!!! Incomprehensible results at isomeric carbon counts:',
                    '\nRegex match: %s\nName processed: %s\n' % (
                        icc, name
                    ),
                    '\n\n\n'
                ]))
        
        return icc
    
    def headgroup_from_lipid_name(self, name, database = 'SwissLipids'):
        """
        For one database record attempts to identify the lipid class
        by looking up keywords.
        Calls greek name identification, greek fatty acid names are
        identified as 'FA'.
        """
        
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
    
    def process(self, name, database = 'swisslipids'):
        """
        The main method of this class. Processes a lipid name string
        and returns a standard name, prefix, carbon counts and
        unsaturations.
        """
        
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
            ccexp = 1 if hg in {'FA', 'MAG'} or lyso else 3 if hg == 'TAG' else 2
            
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


class MoleculeDatabaseAggregator(object):
    
    def __init__(
            self,
            resources = {
                'SwissLipids': (SwissLipids, {}),
                'LipidMaps': (LipidMaps, {})
            },
            tolerance = 20,
            fa_args = None,
            sph_args = None
        ):
        """
        
        """
        
        self.resources = resources
        self.dbs = {}
        self.tolerance = tolerance
        
        self.fa_args  = fa_args or {'c': (4, 36), 'u': (0, 10)}
        self.sph_args = sph_args or {'c': (16, 22), 'u': (0, 1)}
        #self.db = np.array()
    
    def reload(self, children = False):
        
        modname = self.__class__.__module__
        mod = __import__(modname, fromlist=[modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)
    
    def auto_metabolites(
            self,
            fa_args = None,
            sph_args = None,
            sum_only = True
        ):
        
        fa_args  = fa_args  or self.fa_args
        sph_args = sph_args or self.sph_args
        
        masses = []
        data   = []
        
        for name in lipid.sphingolipids:
            
            sys.stdout.write('\t:: Generating `%s`\n' % name)
            
            cls = getattr(lipid, name)
            gen = cls(
                fa_args  = copy.copy(fa_args),
                sph_args = copy.copy(sph_args),
                sum_only = sum_only
            )
            
            for m, d in gen.iterlines():
                
                masses.append(m)
                data.append(d)
        
        self.masses = np.array(masses, dtype = np.float)
        self.data = np.array(data, dtype = np.object)
        
        self.sort()
    
    def sort(self):
        
        self.data = self.data[self.masses.argsort(),:]
        self.masses.sort()
    
    def load(self):
        """
        Loads all databases and generates main array.
        """
        
        for name, (cls, resargs) in self.resources.items():
            
            res = cls(**resargs)
            
            self.dbs[name] = np.array(
                list(res.iterlines())
            )
    
    def ilookup(self, m):
        
        return lookup.findall(self.masses, m, t = self.tolerance)
    
    def lookup(self, m):
        
        i = self.ilookup(m)
        
        return (
            self.masses[i],
            self.data[i]
        )
    
    def lookup_accuracy(self, m):
        
        r = self.lookup(m)
        
        a = np.array([
            (m - rm) / m * 10**6 for rm in r[0]
        ])
        
        return r[0], r[1], a
    
    def adduct_lookup(self, mz, adducts = None, ionm = None, charge = 1):
        
        result = {}
        
        mz = mzmod.Mz(mz)
        
        if ionm in {'pos', 'neg'}:
            
            adducts = list(common.ad2ex[charge][ionm].keys())
        
        methods = dict((ad, common.exact_method[ad]) for ad in adducts)
        
        for ad, method in iteritems(methods):
            
            exmz = getattr(mz, method)()
            
            result[ad] = self.lookup_accuracy(exmz)
        
        return result


class TestMoldb(object):
    
    def __init__(self):
        
        pass
