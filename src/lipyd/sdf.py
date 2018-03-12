#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `lipyd` python module
#
#  Copyright (c) 2014-2018 - EMBL
#
#  File author(s): Dénes Türei (denes@ebi.ac.uk)
#
#  Distributed under the GNU GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://www.ebi.ac.uk/~denes
#

import imp
import sys
import re
import pybel

resyn = re.compile(
    r'(^[A-Z]{2,})\(([0-9]+:[0-9]+)\(.*\)/([0-9]+:[0-9]+)\(.*\)\)'
)
rehg  = re.compile(r'^([A-Z]{2,4})(\(.*\))')
refa  = re.compile(r'C([0-9]+:[0-9]+)n?-?[-0-9]*$')
refa2 = re.compile(r'([0-9]{1,2}:[0-9])\(?[0-9EZ]*\)?$')
hgsyn = {
    'TG': 'TAG',
    'DG': 'DAG'
}

class SdfReader(object):
    
    annots = {
        'PUBCHEM_CID': 'pubchem',
        'SYNONYMS': 'synonym',
        'INCHI': 'inchi',
        'INCHIKEY': 'inchikey',
        'COMMON_NAME': 'commname',
        'SYSTEMATIC_NAME': 'sysname'
    }
    
    def __init__(self, fp, annots = {}, silent = False):
        """
        :param file fp: An open file pointer to the SDF file.
        """
        self.fp = fp
        self.name = self.fp.name
        self.data = {}
        self.mainkey  = {}
        self.indexed = False
        self.silent = silent
        
        self.annots.update(annots)
        
        for annot in self.annots.values():
            
            setattr(self, annot, {})
        
        self._byte_mode()
        self.index()
    
    def reload(self):
        
        modname = self.__class__.__module__
        mod = __import__(modname, fromlist=[modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)
    
    def _byte_mode(self):
        
        if hasattr(self.fp, 'mode'):
            
            if 'b' not in self.fp.mode:
                
                self.fp.close()
                self.fp = open(self.name, 'rb')
    
    def read(self,
            index_only = True,
            one_record = False,
            go_to = 0):
        
        self.fp.seek(-1, 2)
        eof = self.fp.tell()
        self.fp.seek(go_to)
        
        expect_new = True
        molpart = None
        annotpart = False
        annot_or_id = False
        _id = None
        mol = ''
        this_offset = None
        offset = 0
        annot  = {}
        annotkey = None
        
        for l in self.fp:
            
            #print(l)
            
            llen = len(l)
            l = l.decode('utf-8')
            sl = l.strip()
            
            if not molpart:
                
                if annot_or_id and len(sl) and sl[0] != '>' and sl[0] != '$':
                    expect_new = True
                
                if annotkey and annotkey in self.annots:
                    annot[annotkey] = sl
                    annotkey = None
                
                if sl[:3] == '> <':
                    annot_or_id = False
                    annotpart = True
                    annotkey = sl[3:-1]
                
                if annotpart and sl == '':
                    annot_or_id = True
                
                if expect_new and len(l):
                    
                    _id = sl
                    annot = {}
                    this_offset = offset
                    expect_new = False
                    molpart = 1
                    comment = ''
            
            elif molpart == 1:
                
                source = sl
                molpart += 1
            
            elif molpart == 2 and len(sl):
                
                if not sl[0].isdigit():
                    
                    comment = '%s %s' % (comment, sl)
                    
                else:
                    
                    molpart = 3
            
            if sl == '$$$$':
                expect_new = True
            
            if molpart == 3:
                
                if not index_only:
                    
                    mol = '%s%s' % (mol, l)
                
                if sl == 'M  END':
                    molpart = None
                    annotpart = True
                    annot_or_id = True
            
            if expect_new or self.fp.tell() == eof:
                
                if one_record:
                    
                    return {
                        'id': _id,
                        'source': source,
                        'comment': comment,
                        'mol': mol,
                        'annot': annot
                    }
                
                # this is indexing: we build dicts of names
                self.mainkey[_id] = this_offset
                
                if 'COMMON_NAME' in annot:
                    
                    m = refa2.match(annot['COMMON_NAME'])
                    
                    if m:
                        
                        if 'SYNONYMS' not in annot:
                            
                            annot['SYNONYMS'] = 'FA(%s)' % m.groups()[0]
                            
                        else:
                            
                            annot['SYNONYMS'] = '%s;FA(%s)' % (
                                annot['SYNONYMS'],
                                m.groups()[0]
                            )
                
                for k, v in self.annots.items():
                    
                    if k in annot:
                        
                        if k == 'SYNONYMS':
                            
                            syns = set(syn.strip() for syn in annot[k].split(';'))
                            
                            syns2 = set([])
                            
                            for syn in syns:
                                
                                m = rehg.match(syn)
                                
                                if m:
                                    
                                    m = m.groups()
                                    
                                    if m[0] in hgsyn:
                                        
                                        syns2.add('%s%s' % (hgsyn[m[0]], m[1]))
                            
                            syns.update(syns2)
                            syn2 = set([])
                            
                            for syn in syns:
                                
                                m = resyn.match(syn)
                                
                                if m:
                                    
                                    syns2.add('%s(%s/%s)' % m.groups())
                                
                                m = refa.match(syn)
                                
                                if m:
                                    
                                    syns2.add('FA(%s)' % m.groups()[0])
                            
                            syns.update(syns2)
                            
                            for syn in syns:
                                
                                if syn not in self.synonym:
                                    self.synonym[syn] = set([])
                                
                                self.synonym[syn].add(this_offset)
                            
                        else:
                            
                            getattr(self, v)[annot[k]] = this_offset
                
                if not index_only:
                    
                    self.data[this_offset] = {
                        'id': _id,
                        'source': source,
                        'comment': comment,
                        'mol': mol,
                        'annot': annot
                    }
            
            offset += llen
        
        if index_only:
            
            self.indexed = True
    
    def index(self):
        
        self.read(index_only = True)
        self.index_info()
    
    def get_record(self, name, typ):
        
        result = []
        
        if hasattr(self, typ):
            
            index = getattr(self, typ)
            
            if name in index:
                
                if typ == 'synonym':
                    
                    for offset in index[name]:
                        
                        result.append(
                            self.read(
                                index_only = False,
                                one_record = True,
                                go_to = offset
                            )
                        )
                    
                else:
                    
                    offset = index[name]
                    result.append(
                        self.read(
                            index_only = False,
                            one_record = True,
                            go_to = offset
                        )
                    )
        
        return result
    
    def get_obmol(self, name, typ, use_mol = False):
        
        rec = self.get_record(name, typ)
        
        for r in rec:
            
            if use_mol:
                yield self.record_to_obmol_mol(r)
            else:
                yield self.record_to_obmol(r)
    
    def record_to_obmol(self, record):
        
        if 'INCHI' in record['annot']:
            
            return pybel.readstring('inchi', record['annot']['INCHI'])
        
        else:
            
            sys.stdout.write(
                'No InChI for `%s`!\n' % record['annot']['COMMON_NAME']
            )
    
    def record_to_obmol_mol(self, record):
        
        return pybel.readstring('mol', self.get_mol(record))
    
    def get_mol(self, record):
        """
        Returns structure as a string in mol format.
        """
        
        return '%s\n  %s\n%s\n%s' % (
            record['id'],
            record['source'],
            record['comment'],
            record['mol']
        )
    
    def write_mol(self, name, typ, outf = None, return_data = False):
        
        outf = outf or '%s_%s_%s.mol'
        
        rr = self.get_record(name, typ)
        
        if not rr:
            
            return None
        
        if type(rr) is not list:
            
            rr = [rr]
        
        for r in rr:
            
            _outf = outf % (
                name.replace('/', '.'),
                r['annot']['COMMON_NAME'].replace('/', '.').replace(' ', '..')
                    if 'COMMON_NAME' in r['annot']
                    else '',
                r['id']
            )
            
            r['molfile'] = _outf
            
            with open(_outf, 'w') as fp:
                
                _ = fp.write(
                    self.get_mol(r)
                )
        
        if return_data:
            
            return rr
    
    def index_info(self):
        
        if not self.silent:
            
            sys.stdout.write('\t:: Indexed %u records from `%s`.\n' % (
                len(self.mainkey),
                self.name
            ))
    
    def __del__(self):
        
        self.fp.close()
