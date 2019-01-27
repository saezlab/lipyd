#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `lipyd` python module
#
#  Copyright (c) 2014-2019 - EMBL
#
#  File author(s):
#  Dénes Türei (turei.denes@gmail.com)
#  Igor Bulanov
#
#  Distributed under the GNU GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://denes.omnipathdb.org/
#

class Recalibration(object):
    
    
    def __init__(
            self,
            first = None,
            last = None,
            by_sample = None,
            **kwargs,
        ):
        """
        Provides recalibration data for ``lipyd.sample.Sample`` and
        ``lipyd.sample.SampleSet`` objects.
        
        :arg float,numpy.ndarray first:
            Error in ppm for the first sample. Either a single float or an
            array of floats with the same length as number of features in
            the sample or an array of 2 columns: first with m/z values,
            second with errors in ppm above the corresponding m/z value.
        :arg float,numpy.ndarray last:
            Error in ppm for the last sample. Same as ``first``. If not
            provided all samples will be corrected by ``first``.
        :arg dict by_sample:
            Errors for each sample. Keys are sample IDs, values are as
            detailed at argument ``first``: float or array with one or
            two columns. ``**kwargs`` handled the same way as ``by_sample``
            but the latter has priority.
        """
        
        self.first = first
        self.last = last
        self.by_sample = by_sample or kwargs
    
    
    def recalibrate(self, samples):
        """
        Performs correction on m/z values in a ``lipyd.sample.Sample`` or a
        ``lipyd.sample.SampleSet`` object.
        
        :arg Sample,SampleSet samples:
            A ``lipyd.sample.Sample`` or a ``lipyd.sample.SampleSet`` object.
        """
        
        self.samples = samples
        self.numof_samples = self.samples.numof_samples
        self.numof_features = len(self.samples)
    
    
    def ppm_array(self):
        
        pass
