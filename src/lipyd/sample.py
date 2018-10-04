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
import numpy as np


class FeatureAttributes(object):
    
    def __init__(self):
        """
        
        """
        
        pass
    
    def reload(self):
        
        modname = self.__class__.__module__
        mod = __import__(modname, fromlist=[modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)
    
    def __len__(self):
        
        return self.data.shape[0]
    
    def charges(self):
        """
        Returns a set of ion charges observed in the sample.
        """
        
        return sorted(set(self.charge))


class Sample(object):
    
    def __init__(
            self,
            mzs,
            intensities = None,
            rts = None,
            attrs = None,
        ):
        """
        Represents one LC MS/MS run.
        Has at least a vector of m/z's, optionally vector of intensites,
        retention times and other metadata.
        """
        
        self.var     = set()
        self.missing = set()
        
        if mzs is None:
            
            raise ValueError('Sample object must have at least m/z values.')
        
        self.add_var(mzs, 'mzs')
        self.add_var(intensities, 'intensities')
        self.add_var(rts, 'rts')
        
        self.attrs = attrs or {}
    
    def reload(self):
        
        modname = self.__class__.__module__
        mod = __import__(modname, fromlist=[modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)
    
    def add_var(self, data, attr):
        
        if attr != 'mzs' and data is not None and len(data) != len(self):
            
            raise ValueError(
                'Sample object: length of `%s` (%u) is not the same '
                'as number of features in the sample (%u).' % (
                    attr, len(data), len(self),
                )
            )
        
        setattr(self, attr, data)
        
        if data is None:
            
            self.missing.add(attr)
            
        else:
            
            self.var.add(attr)
    
    def sort_all(self, by = 'mzs', desc = False, resort = False):
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
        
        return len(self.mzs)
