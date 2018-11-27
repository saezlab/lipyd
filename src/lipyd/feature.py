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

# hello!

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


class PeakSize(FeatureAnalyzer):
    
    def __init__(
            self,
            samples,
            protein_samples,
            threshold = 2.0,
            name = 'peaksize',
        ):
        
        self.threshold = threshold
        
        FeatureAnalyzer.__init__(
            self,
            name = name,
            samples = samples,
            method = self.peak_size,
            protein_samples = protein_samples,
        )
    
    def peak_size(self, intensities, protein_samples, **kwargs):
        
        protein   = intensities[protein_samples]
        noprotein = intensities[np.logical_not(protein_samples)]
        
        if np.any(np.isnan(protein)):
            
            return False
        
        if np.all(np.isnan(noprotein)):
            
            return True
        
        return np.min(protein) > np.nanmax(noprotein) * self.threshold


class Slope(FeatureAnalyzer):
    
    def __init__(
            self,
            samples,
            protein_profile,
            protein_samples,
            highest_samples,
            name = 'slope',
        ):
        
        self.threshold = threshold
        
        FeatureAnalyzer.__init__(
            self,
            name = name,
            samples = samples,
            method = self.slope,
            protein_samples = protein_samples,
        )
    
    def slope(
            self,
            intensities,
            protein_profile,
            protein_samples,
            highest_samples,
        ):
        
        pass
