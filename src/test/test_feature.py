#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `lipyd` python module
#
#  Copyright (c) 2015-2018 - EMBL
#
#  File author(s):
#  Dénes Türei (turei.denes@gmail.com)
#  Igor Bulanov
#
#  Distributed under the GNU GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://www.ebi.ac.uk/~denes
#

import pytest

import numpy as np

import lipyd.sample as sample
import lipyd.feature as feature
import lipyd.sampleattrs as sampleattrs


class TestFeature(object):
    """ """
    
    @pytest.fixture(autouse = True)
    def auto_inject_fixture(self):
        """ """
        
        sample_ids = ['A10', 'A11', 'A12', 'B1']
        
        intensities = np.vstack([
            [.1] * 10,
            np.arange(10) * .1,
            (np.arange(10) * .1)[::-1],
            [.1] * 10,
        ])
        intensities = intensities.transpose()
        
        mzs = np.random.random(10)
        
        self.samples = sample.SampleSet(
            mzs = mzs,
            intensities = intensities,
            ionmode = 'pos',
            sample_ids = sample_ids,
            sample_id_proc = sampleattrs.plate_sample_id_processor(),
        )
        
        self.protein_samples = self.samples.get_selection(['A11', 'A12'])
    
    def test_feature_analyzer(self):
        """ """
        
        def method(intensities, protein_samples, **kwargs):
            """

            Parameters
            ----------
            intensities :
                
            protein_samples :
                
            **kwargs :
                

            Returns
            -------

            """
            
            return (
                np.nanmin(intensities[protein_samples]) >
                np.nanmax(intensities[np.logical_not(protein_samples)])
            )
        
        fea = feature.FeatureAnalyzer(
            name = 'peaksize',
            samples = self.samples,
            method = method,
            protein_samples = self.protein_samples.selection,
        )
        
        assert np.all(
            self.samples.feattrs.peaksize ==
            np.array([
                False, False,  True,  True,  True,
                True,  True,  True, False, False
            ])
        )
    
    def test_peak_size(self):
        """ """
        
        fea = feature.PeakSize(
            samples = self.samples,
            protein_samples = self.protein_samples.selection,
        )
        
        assert np.all(
            self.samples.feattrs.peaksize ==
            np.array([
                False, False, False,  True,  True,
                True,  True, False, False, False
            ])
        )
