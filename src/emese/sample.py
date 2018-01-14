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
import numpy as np

import emese.mgf as mgf


class Sample(mgf.MgfReader):
    
    def __init__(self, *args, **kwargs):
        
        self.sorted_by = None
        
        mgf.MgfReader.__init__(self, *args, **kwargs)
        
        self.read()
        self.sort_all()
    
    def reload(self):
        modname = self.__class__.__module__
        mod = __import__(modname, fromlist=[modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)
    
    def sort_all(self, by = 'mz', desc = False, resort = False):
        """
        Sorts all data arrays according to values in one of them.
        """
        
        if by == self.sorted_by and not resort:
            
            return None
        
        isort = np.argsort(getattr(self, by))
        
        if desc:
            isort = isort[::-1]
        
        for var in self.var:
            
            setattr(self, var, getattr(self, var)[isort])
        
        self.sorted_by = by
    
    def __len__(self):
        """
        Returns number of MS1 ions (m/z's) detected in the sample.
        """
        
        return len(self.mz)
    
    def charges(self):
        """
        Returns a set of ion charges observed in the sample.
        """
        
        return sorted(set(self.charge))
