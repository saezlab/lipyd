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

import warnings
import numpy as np

import lipyd.sample as sample
import lipyd.sampleattrs as sampleattrs
import lipyd.settings as settings


class TestSample(object):
    """ """
    
    @pytest.fixture(autouse = True)
    def auto_inject_fixture(self):
        """ """
        
        sample_ids = ['A12', 'A11', 'B2', 'B1']
        
        mzs = np.arange(40)
        mzs.shape = (10, 4)
        
        rts = np.arange(40)[::-1]
        rts.shape = (10, 4)
        
        self.samples = sample.SampleSet(
            mzs = mzs,
            rts = rts,
            ionmode = 'pos',
            sample_ids = sample_ids,
            sample_id_proc = sampleattrs.plate_sample_id_processor(),
        )
    
    def test_feature_idx(self):
        """ """
        
        idx = sample.FeatureIdx(10)
        
        assert np.all(idx._current == idx._original)
        
        rand10 = np.empty((10,))
        sort10 = rand10.argsort()
        idx.sort_all(sort10)
        
        assert np.all(idx._current == sort10)
        # the smallest number is at index zero currently
        # as arrays are sorted by the data
        assert idx._original[rand10.argmin()] == 0
        
        idx.sort_all()
        
        assert np.all(idx._current == idx._original)
    
    def test_feature_base(self):
        """ """
        
        a0 = np.random.random(10)
        a1 = np.random.random(10)
        a2 = np.random.random(10)
        b1 = np.random.random(10)
        
        a0min = a0.argmin()
        a2a0min = a2[a0min]
        b1a0min = b1[a0min]
        
        f0 = sample.FeatureBase(a = a0)
        f1 = sample.FeatureBase(a = a1, b = b1, sorter = f0.sorter)
        f2 = sample.FeatureBase(a = a2, sorter = f0.sorter)
        
        f0.sort_all('a')
        
        assert f2.a[0] == a2a0min
        assert f1.b[0] == b1a0min
    
    def test_feature_base_3d(self):
        """ """
        
        a0 = np.random.random((10, 7, 5))
        a1 = np.random.random((10, 7, 5))
        a2 = np.random.random((10, 7, 5))
        b1 = np.random.random((10, 7, 5))
        
        a0min = a0[:,0,0].argmin()
        a2a0min = a2[a0min,0,0]
        b1a0min = b1[a0min,0,0]
        
        f0 = sample.FeatureBase(a = a0)
        f1 = sample.FeatureBase(a = a1, b = b1, sorter = f0.sorter)
        f2 = sample.FeatureBase(a = a2, sorter = f0.sorter)
        
        f0.sort_all('a', indices = (0, 0))
        
        assert f2.a[0,0,0] == a2a0min
        assert f1.b[0,0,0] == b1a0min
        
        with pytest.warns(UserWarning):
            
            f2.sort_all('a')
    
    def test_sampleset(self):
        """ """
        
        assert self.samples.attrs.attrs[0].sample_id == ('A', 12)
    
    def test_samplesorter(self):
        """ """
        
        samples = self.samples
        
        # repeating these as pytest makes a copy of the object
        order = ['A11', 'A12', 'B1', 'B2']
        samples.sort_by_sample_ids(order)
        
        assert samples.mzs_by_sample[0, 1] == 0
        assert samples.rts[4, 3] == 21
        
        samples.sort_all(by = 'mzs', desc = True)
        
        assert samples.mzs_by_sample[0, 1] == 36
        assert samples.rts[4, 3] == 17
    
    def test_sampledata(self):
        """ """
        
        samples = self.samples
        
        # repeating these as pytest makes a copy of the object
        order = ['A11', 'A12', 'B1', 'B2']
        samples.sort_by_sample_ids(order)
        samples.sort_all(by = 'mzs', desc = True)
        
        data0 = ['a',  'b',  'c',   'd']
        order = ['B2', 'B1', 'A12', 'A11']
        
        sdata = sampleattrs.SampleData(
            var0 = data0,
            sample_ids = order,
            samples = samples,
            sample_id_proc = sampleattrs.plate_sample_id_processor(),
        )
        
        assert sdata.var0[0] == 'd' and sdata.var0[-1] == 'a'
        assert samples.mzs_by_sample[9,1] == 0
        assert (
            samples.attrs.sample_index_to_id ==
            sdata.attrs.sample_index_to_id
        )
        assert samples.attrs.sample_index_to_id[2] == ('B', 1)
        assert samples.rts[0,2] == 0
        assert id(sdata) in samples._sample_data
        assert id(samples) in sdata._sample_data
    
    def test_sampleset_from_peaks(self):
        """ """
        
        peaksfile = settings.get('peaks_example')
        
        reader = sample.SampleReader(
            input_type = 'peaks',
            fname = peaksfile
        )
        
        samples = reader.get_sampleset(
            sampleset_args = {
                'sample_id_proc': sampleattrs.plate_sample_id_processor(),
            }
        )
        
        assert abs(samples.mzs_by_sample[7,3] - 375.0018) < 0.0001
        assert samples.attrs.sample_index_to_id[-1] == ('A', 12)
        assert samples.attrs.attrs[0].attrs['label']['sample_id'] == ('A', 6)
    
    def test_sampleselection(self):
        """ """
        
        samples = self.samples
        order = ['A11', 'A12', 'B1', 'B2']
        samples.sort_by_sample_ids(order)
        
        sel = sampleattrs.SampleSelection(
            selection = ['A12', 'B1'],
            samples = samples,
        )
        
        assert np.all(sel.selection == np.array([False, True, True, False]))
        
        samples.sort_by_sample_ids(['A11', 'B2', 'A12', 'B1'])
        
        assert np.all(sel.selection == np.array([False, False, True, True]))
        
        sel = sampleattrs.SampleSelection(
            selection = np.array([True, True, False, False]),
            samples = samples
        )
        
        assert np.all(sel.selection == np.array([True, True, False, False]))
        
        samples.sort_by_sample_ids(['A12', 'B1', 'A11', 'B2'])
        
        assert np.all(sel.selection == np.array([False, False, True, True]))
    
    def test_get_selection(self):
        """ """
        
        samples = self.samples
        order = ['A11', 'A12', 'B1', 'B2']
        samples.sort_by_sample_ids(order)
        
        # SampleSelection
        
        sel = samples.get_selection(selection = ['A12', 'B1'])
        
        assert np.all(sel.selection == np.array([False, True, True, False]))
        
        samples.sort_by_sample_ids(['A11', 'B2', 'A12', 'B1'])
        
        assert np.all(sel.selection == np.array([False, False, True, True]))
        
        sel = samples.get_selection(np.array([True, True, False, False]))
        
        assert np.all(sel.selection == np.array([True, True, False, False]))
        
        samples.sort_by_sample_ids(['A12', 'B1', 'A11', 'B2'])
        
        assert np.all(sel.selection == np.array([False, False, True, True]))
        
    def test_get_sample_data(self):
        """ """
        
        samples = self.samples
        order = ['A11', 'A12', 'B1', 'B2']
        samples.sort_by_sample_ids(order)
        
        # SampleData
        
        dt = samples.get_sample_data(data0 = np.array([7, 77, 777, 7777]))
        
        assert np.all(dt.data0 == np.array([7, 77, 777, 7777]))
        
        samples.sort_by_sample_ids(['A11', 'B2', 'A12', 'B1'])
        
        assert np.all(dt.data0 == np.array([7, 7777, 77, 777]))
