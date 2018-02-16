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

from future.utils import iteritems

import os
import sys
import imp
import re
import struct
import itertools
import collections

import numpy as np

try:
    import pybel
except:
    sys.stdout.write('\t:: No module `pybel` available.\n')

import lipyd._curl as _curl
import lipyd.common as common
import lipyd.settings as settings
import lipyd.mz as mz
import lipyd.progress as progress
import lipyd.sdf as sdf


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
        
        with open(self.fname, 'w') as fpe:
            
            while True:
                
                block = fpa.read(4096)
                if not block:
                    break
                
                fpe.write(block.decode('utf-8'))
        
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
    
    def __init__(self, levels = set(['Species']), silent = False):
        
        self.silent = silent
        self.levels = levels
        self.url = settings.get('swisslipids_url')
        self.load()
        self.make_index()
    
    def load(self):
        
        self.close_gzfile()
        
        self._curl = _curl.Curl(self.url, silent = False,
                                compr = 'gz', large = True)
        self._gzfile = self._curl.result
    
    def iterfile(self):
        
        self._gzfile.fileobj.seek(0)
        _ = self._gzfile.readline()
        
        for line in self._gzfile:
            
            yield line
    
    @staticmethod
    def names(line):
        
        return '|'.join(line[2:5])
    
    def make_index(self):
        
        self.close_plainfile()
        
        self.load()
        self.index = collections.defaultdict(lambda: set([]))
        
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
                
                ll = l.decode('utf-8').split('\t')
                
                if ll[1] in self.levels:
                    
                    names = self.names(ll)
                    self.index[ll[0]].add(offset) # SwissLipids ID
                    self.index[ll[8]].add(offset) # SMILES
                    self.index[ll[10]].add(offset) # InChI key
                    
                    for n in names.split('|'):
                        
                        self.index[n].add(offset)
                
                offset = self._gzfile.tell()
                fpp.write(l)
        
        if not self.silent:
            self.prg.terminate()
        
        self.index = dict(self.index)
        self._plainfile = open(self._plainfilename, 'r')
    
    def get_record(self, name):
        
        result = []
        
        if name in self.index:
            
            for offset in self.index[name]:
                
                self._plainfile.seek(offset)
                result.append(self._plainfile.readline().strip().split('\t'))
        
        return result
    
    def get_obmol(self, name):
        
        return [self.to_obmol(rec) for rec in self.get_record(name)]
    
    @staticmethod
    def to_obmol(record):
        
        if not record[8]:
            
            return None
        
        mol = pybel.readstring('smi', record[8])
        mol.title = '|'.join(record[2:5])
        mol.chebi = record[24] if len(record) > 24 else ''
        mol.lipidmaps = record[25] if len(record) > 25 else ''
        mol.hmdb = record[26] if len(record) > 26 else ''
        mol.smiles = record[8]
        mol.swl_exact_mass = float(record[14]) if record[14] else None
        mol.swl_formula = record[11]
        mol.inchikey = (
            record[10].split('=')[1]
            if '=' in record[10]
            else record[10])
        mol.level = record[1]
        
        return mol
    
    def __iter__(self):
        
        if not self.silent:
            
            self._gzfile.fileobj.seek(-4, 2)
            ucsize = struct.unpack('I', self._gzfile.fileobj.read(4))[0]
            self.prg = progress.Progress(ucsize, 'Processing SwissLipids', 101)
            self._gzfile.fileobj.seek(0)
        
        nosmiles = 0
        
        for line in self.iterfile():
            
            if not self.silent:
                
                self.prg.step(len(line))
            
            line = line.decode('utf-8').strip().split('\t')
            
            if len(line) > 22 and line[1] in self.levels:
                
                if not line[8]:
                    nosmiles += 1
                    continue
                
                self.to_obmol(line)
                
                yield mol
        
        if not self.silent:
            
            self.prg.terminate()
        
        self._curl.close()
    
    def __del__(self):
        
        self.close_gzfile()
        self.close_plainfile()
    
    def close_gzfile(self):
        
        if hasattr(self, '_gzfile'):
            
            self._gzfile.close()
    
    def close_plainfile(self):
        
        if hasattr(self, '_plainfile'):
            
            self._plainfile.close()


class LipidNameProcessor(object):
    
    def __init__(self):
        
        self.lipnamesf = settings.get('lipnamesf')
        self.adducts_constraints = settings.get('adducts_constraints')
        self.recount1 = re.compile(r'\(([Odt]?)-?([0-9]{1,2}):([0-9]{1,2})\)')
        self.recount2 = re.compile(r'\(([Odt]?)-?([0-9]{1,2}):([0-9]{1,2})/'
                                   r'([Odt]?)-?([0-9]{1,2}):([0-9]{1,2})/?'
                                   r'([Odt]?)-?([0-9]{0,2}):?([0-9]{0,2})\)')
        self.recount3 = re.compile(r'\(([Odt]?)-?([0-9]{1,2}):([0-9]{1,2})/?'
                                   r'([Odt]?)-?([0-9]{0,2}):?([0-9]{0,2})/?'
                                   r'([Odt]?)-?([0-9]{0,2}):?([0-9]{0,2})\)')
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
    
    def carbon_counts(self, name):
        
        # regex finds the total carbon count
        cc1 = self.recount1.findall(name)
        # regex finds 2-3 fatty acids
        cc2 = self.recount3.findall(name)
        
        # the total carbon count
        ccpart = (
            [cc1[0][0], int(cc1[0][1]), int(cc1[0][2])]
            if len(cc1) else
            ['', np.nan, np.nan]
        )
        
        # carbon counts of fatty acids
        if len(cc2s) and (
            any(map(lambda cc2: cc2[4], cc2s)) or
            swl_parsed == 'FA' or
            lyso
        ):
            
            faccparts = (
                list(
                    map(
                        lambda cc2:
                            [
                                # FA1
                                cc2[0],
                                int(cc2[1]),
                                int(cc2[2]),
                                # FA2
                                cc2[3],
                                int(cc2[4]) if cc2[4] else np.nan,
                                int(cc2[5]) if cc2[5] else np.nan,
                                # FA3
                                cc2[6],
                                int(cc2[7]) if cc2[7] else np.nan,
                                int(cc2[8]) if cc2[8] else np.nan
                            ],
                        filter(
                            lambda cc2:
                                # if this is a Lyso species
                                # or single fatty acid
                                # we have only one cc:unsat
                                # otherwise we must have at least 2
                                cc2[4] or swl_parsed == 'FA' or lyso,
                            cc2s
                        )
                    )
                )
            )
            
        else:
            faccparts = [
                [
                    '', np.nan, np.nan,
                    '', np.nan, np.nan,
                    '', np.nan, np.nan
                ]
            ]
        
        for faccpart in faccparts:
            
            if cc2s and not cc1:
                
                ccpart = [
                    faccpart[0] or faccpart[3] or faccpart[6],
                    np.nansum([faccpart[1], faccpart[4], faccpart[7]]),
                    np.nansum([faccpart[2], faccpart[5], faccpart[8]])
                ]
            
            # making sure ethers are idenfified:
            if 'O' not in ccpart[0] and ('O-' in lip or '-O' in lip):
                
                ccpart[0] = 'O%s' % ccpart[0]
            
            if cc1 and 'O' not in cc1[0][0] and 'O' in ccpart[0]:
                
                cc1[0] = ('O%s' % cc1[0][0], cc1[0][1], cc1[0][2])
        
        return cc1, cc2
    
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
                print(shortname, matched)
                if sum(matched) == len(kwset['pos']) and \
                    sum(matched) > 0:
                    matched = [kw in name for kw in kwset['neg']]
                    if sum(matched) == 0:
                        ether = '(O-' in name
                        return shortname, spec['pos_adduct'], spec['neg_adduct'], ether
        
        return self.process_fa_name(name)
    
    def process_fa_name(self, name):
        """
        Identifies fatty acids based on their greek name.
        """
        
        if name in self.fa_greek:
            return 'FA', None, None, False
        
        return None, None, None, None
    
    def gen_fa_greek(self):
        
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
        
        self.fa_greek = {}
        
        for cc, uns, end in itertools.product(
            fa_greek_parts['cc'].items(),
            fa_greek_parts['uns'].items(),
            fa_greek_parts['end'].items()):
            
            if len(uns[0]) and end[1] == 0:
                continue
            
            self.fa_greek['%s%s%s' % (cc[0], uns[0], end[0])] = (cc[1], uns[1] * end[1])
