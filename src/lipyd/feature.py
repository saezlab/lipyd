#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `lipyd` python module
#
#  Copyright (c) 2015-2018 - EMBL
#
#  File author(s): Dénes Türei (turei.denes@gmail.com)
#
#  Distributed under the GNU GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://www.ebi.ac.uk/~denes
#

from future.utils import iteritems
from past.builtins import xrange, range

import numpy as np

import lipyd.sampleattrs as sampleattrs


class FeatureAnalyzer(object):
    
    def __init__(self, name, samples, method, **variables):
        
        self.name = name
        self.samples = samples
        self.method = method
        self.variables = variables
        
        result = self.run()
        self.samples.feattrs._add_var(result, self.name)
    
    def run(self):
        """
        Applies the method for each feature and returns an array of the
        results.
        """
        
        result = []
        
        for i in xrange(len(self.samples)):
            
            featurevars = {}
            
            for var in self.samples.feattrs.var:
                
                featurevars[var] = getattr(self.samples.feattrs, var)[i]
            
            for var in self.samples.var:
                
                # safe to do this as first axis is always the features
                featurevars[var] = getattr(self.samples, var)[i]
            
            result.append(self.method(**featurevars, **self.variables))
        
        return np.array(result)


class ProfileFeatureAnalyzer(FeatureAnalyzer):
    
    def __init__(self, samples):
        
        FeatureAnalyzer.__init__(self, samples = samples)
