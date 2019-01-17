#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `lipyd` python module
#
#  Copyright (c) 2015-2019 - EMBL
#
#  File author(s): Dénes Türei (turei.denes@gmail.com)
#
#  Distributed under the GNU GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://denes.omnipathdb.org/
#


import os
import imp
import re
import regex
import itertools


class MfqlReader(object):
    """ """
    
    section_names = {
        'QUERYNAME',
        'DEFINE',
        'IDENTIFY',
        'SUCHTHAT',
        'REPORT',
    }
    
    precedence = (
        re.compile(r'(WITH|WHERE)(?![^(]*\))'),
        re.compile(r'(AND|OR)(?![^(]*\))'),
        re.compile(r'(=+)(?![^(]*\))'),
        re.compile(r'[\s]([-\+])(?![^(]*\))'),
        re.compile(r'(,)(?![^(]*\))'),
        re.compile(r'[^"](%)[^\s](?![^(]*\))'),
        re.compile(r'\s(IN|in)\s(?![^(]*\))')
    )
    
    retitle = re.compile(r'(%s)' % '|'.join(section_names))
    reblank = re.compile(r'\s+')
    renewln = re.compile(r'[\n\r]+')
    repar   = re.compile(r'\((.*)\)')
    recomm  = re.compile(r'^[\s]*(?:#.*)?[\s]*$', flags = 8)
    rerepar = regex.compile(
        r'(?:\(((?>[^\(\)]+|(?R))*)\))?'
        r'([^\(^\)]*)'
    )
    refun   = re.compile(r'(\w+)\((.*)\)')
    
    def __init__(self, fname = None):
        
        self.fname = fname
    
    def reload(self):
        """ """
        
        modname = self.__class__.__module__
        mod = __import__(modname, fromlist=[modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)
    
    def read_raw(self):
        """ """
        
        sections = {}
        
        with open(self.fname, 'r') as fp:
            
            raw_sections = self.retitle.split(fp.read())
            
            for title, content in zip(
                itertools.islice(raw_sections, 1, None, 2),
                itertools.islice(raw_sections, 2, None, 2),
            ):
                
                # removing comment lines
                content = self.recomm.sub(' ', content)
                # removing all new lines and indentation from content
                content = self.reblank.sub(' ', content)
                
                tokens = self.tokenize(content)
                
                sections[title] = tokens
        
        self.sections = sections
    
    def read_sections(self):
        """ """
        
        result = {}
        
        for line in self.read_raw():
            
            _title = line.split(maxsplit = 1)[0].split('=', maxsplit = 1)[0]
    
    @classmethod
    def tokenize(cls, s):
        """

        Parameters
        ----------
        s :
            

        Returns
        -------

        """
        
        s = s.strip().strip(';')
        
        if ';' in s:
            
            return tuple(cls.tokenize(si) for si in s.split(';'))
        
        if s and '(' == s[0]:
            
            par = cls.rerepar.findall(s)
            
            return tuple(
                (
                    ('PAR', cls.tokenize(inpar)),
                    cls.tokenize(postpar)
                )
                for inpar, postpar
                in par[:-1]
            )
        
        for sep in cls.precedence:
            
            splt = sep.split(s)
            
            if len(splt) < 3:
                
                continue
            
            for elem0, rel, elem1 in zip(
                itertools.islice(splt, 0, None, 3),
                itertools.islice(splt, 1, None, 3),
                itertools.islice(splt, 2, None, 3),
            ):
                
                return (rel.upper(), cls.tokenize(elem0), cls.tokenize(elem1))
        
        func = cls.refun.match(s)
        
        if func:
            
            func = func.groups()
            
            print(func)
            
            return (func[0].upper(), cls.tokenize(func[1]))
        
        # no more processing necessary
        return s
    
    def queryname(self, line):
        """

        Parameters
        ----------
        line :
            

        Returns
        -------

        """
        
        pass
