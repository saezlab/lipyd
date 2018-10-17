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

class SampleSelection(object):
    
    def __init__(self):
        
        pass


class FeatureAnalyzer(object):
    
    def __init__(self, samples):
        
        self.samples = samples


class ProfileFeatureAnalyzer(FeatureAnalyzer):
    
    def __init__(self, samples):
        
        FeatureAnalyzer.__init__(self, samples = samples)
