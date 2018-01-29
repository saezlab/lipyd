#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `emese` python module
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

from future.utils import iteritems

import os
import sys
import imp
import struct

try:
    import pybel
except:
    sys.stdout.write('\t:: No module `pybel` available.\n')

import emese._curl as _curl
import emese.common as common
import emese.settings as settings
import emese.mz as mz
import emese.progress as progress
import emese.sdf as sdf


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
        self.url = common.get_param('swisslipids_url')
    
    def load(self):
        
        self._curl = _curl.Curl(self.url, silent = False,
                                compr = 'gz', large = True)
        self._gzfile = self._curl.result
        
        if not self.silent:
            
            self._gzfile.fileobj.seek(-4, 2)
            ucsize = struct.unpack('I', self._gzfile.fileobj.read(4))[0]
            self.prg = progress.Progress(ucsize, 'Processing SwissLipids', 101)
            self._gzfile.fileobj.seek(0)
    
    def iterfile(self):
        
        self._gzfile.fileobj.seek(0)
        _ = self._gzfile.readline()
        
        for line in self._gzfile:
            
            yield line
    
    @staticmethod
    def names(line):
        
        return '|'.join(line[2:5])
    
    def __iter__(self):
        
        nosmiles = 0
        
        for line in self.iterfile():
            
            if not self.silent:
                
                self.prg.step(len(line))
            
            line = line.decode('utf-8').strip().split('\t')
            
            if len(line) > 22 and line[1] in self.levels:
                
                if not line[8]:
                    nosmiles += 1
                    continue
                
                mol = pybel.readstring('smi', line[8])
                mol.title = self.names(line)
                mol.chebi = line[24] if len(line) > 24 else ''
                mol.lipidmaps = line[25] if len(line) > 25 else ''
                mol.hmdb = line[26] if len(line) > 26 else ''
                mol.smiles = line[8]
                mol.swl_exact_mass = float(line[14]) if line[14] else None
                mol.swl_formula = line[11]
                mol.inchikey = (
                    line[10].split('=')[1]
                    if '=' in line[10]
                    else line[10])
                mol.level = line[1]
                
                yield mol
        
        if not self.silent:
            
            self.prg.terminate()
        
        self._curl.close()
