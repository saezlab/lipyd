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

from past.builtins import xrange, range

import imp
import re
import collections
import mimetypes

import numpy as np

import lipyd.settings as settings
import lipyd.common as common
import lipyd.reader.xls as xls


refrac = re.compile(r'([A-Z])([0-9]{1,2})')


Fraction = collections.namedtuple('Fraction', ['label', 'start', 'end'])


class SECProfile(object):
    
    def __init__(self, path):
        
        self.path = path
        self.read()
    
    def reload(self, children = False):
        
        modname = self.__class__.__module__
        mod = __import__(modname, fromlist=[modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)
    
    def read(self):
        
        self.guess_format()
        
        if self.format == 'asc':
            
            self.read_asc()
            
        elif self.format == 'xls':
            
            self.read_xls()
    
    def read_asc(self):
        """
        Reads SEC UV absorbance profile from asc file output produced by
        the Unicorn software from GE Healthcare.
        """
        
        start      = None
        end        = None
        frac       = None
        volume     = []
        absorbance = []
        fractions  = []
        
        with open(self.path, 'r') as fp:
            
            for l in fp:
                
                l = l.strip().split('\t')
                
                if len(l) < 2:
                    
                    continue
                
                vol = common.to_float(l[0])
                ab_ = common.to_float(l[1])
                
                if vol and ab_:
                    
                    volume.append(vol)
                    absorbance.append(ab_)
                
                if len(l) > 3:
                    
                    start = end
                    end   = common.to_float(l[2])
                    
                    if start and end and frac:
                        
                        fractions.append(Fraction(frac, start, end))
                    
                    m = refrac.search(l[3])
                    
                    frac = m.groups() if m else None
        
        self.volume = np.array(volume)
        self.absorbance = np.array(absorbance)
        self.fractions = fractions
    
    def read_xls(self):
        """
        Reads SEC UV absorbance profile from MS Excel XLS file output
        produced by ???.
        """
        
        volume     = []
        absorbance = []
        
        tab = xls.read_xls(self.path)
        
        for l in tab:
            
            vol = common.to_float(l[0])
            ab_ = common.to_float(l[1])
            
            if vol and ab_:
                
                volume.append(vol)
                absorbance.append(ab_)
        
        self.volume = np.array(volume)
        self.absorbance = np.array(absorbance)
    
    def auto_fractions(
            self,
            start_volume = .615,
            size = .15,
            start_row = 'A',
            start_col = 5,
            length = 9,
        ):
        
        fractions = []
        
        for i in xrange(length):
            
            start = start_volume + size * i
            end   = start + size
            well  = (ord(start_row) - 65) * 12 + start_col + i - 1
            row   = chr(well // 12 + 65)
            col   = well % 12 + 1
            
            fractions.append(Fraction((row, col), start, end))
        
        return fractions
    
    def guess_format(self):
        
        mime = mimetypes.guess_type(self.path)
        self.format = (
            'xls' if 'excel' in mime or 'openxml' in mime else 'asc'
        )
