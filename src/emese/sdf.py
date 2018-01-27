#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `emese` python module
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

class SdfReader(object):
    
    annots = {
        'PUBCHEM_CID': 'pubchem',
        'SYNONYMS': 'synonym',
        'INCHIKEY': 'inchikey',
        'COMMON_NAME': 'commname',
        'SYSTEMATIC_NAME': 'sysname'
    }
    
    def __init__(self, fp, annots = {}, silent = False):
        """
        :param file fp: An open file pointer to the SDF file.
        """
        self.fp = fp
        self.name = fp.name
        self.data = {}
        self.mainkey  = {}
        self.indexed = False
        
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
            
            l = l.decode('ascii')
            sl = l.strip()
            
            if not molpart:
                
                if annot_or_id and len(sl) and sl[0] != '>' and sl[0] != '$':
                    expect_new = True
                
                if annotkey and annotkey in self.annots:
                    annot[annotkey] = sl
                    annotkey = None
                
                if sl[:3] == '> <':
                    if not annotpart:
                        annot = {}
                    annot_or_id = False
                    annotpart = True
                    annotkey = sl[3:-1]
                
                if annotpart and sl == '':
                    annot_or_id = True
                
                if expect_new and len(l):
                    
                    _id = sl
                    print(_id)
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
                
                for k, v in self.annots.items():
                    
                    if k in annot:
                            
                        if k == 'SYNONYMS':
                            
                            for syn in annot[k].split(';'):
                                
                                self.synonym[syn.strip()] = this_offset
                            
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
            
            offset += len(l)
        
        if index_only:
            
            self.indexed = True
    
    def index(self):
        
        self.read(index_only = True)
        self.index_info()
    
    def get_record(self, name, typ):
        
        if hasattr(self, typ):
            
            index = getattr(self, typ)
            
            if name in index:
                
                offset = index[name]
                
                return self.read(index_only = False, one_record = True, go_to = offset)
    
    def write_mol(self, name, typ, outf = None):
        
        outf = outf or '%s.mol' % name
        
        r = self.get_record(name, typ)
        
        with open(outf, 'w') as fp:
            
            fp.write(
                '%s\n  %s\n%s\n%s' % (
                    r['id'], r['source'], r['comment'], r['mol']
                )
            )
    
    def index_info(self):
        
        sys.stdout.write('\t:: Indexed %u records from `%s`.\n' % (
            len(self.mainkey),
            self.name
        ))
    
    def __del__(self):
        
        self.fp.close()
