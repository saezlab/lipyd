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

import pytest

import re
import numpy as np

import lipyd.sampleattrs as sampleattrs
import lipyd.sample as sample
import lipyd.settings as settings
from lipyd.common import basestring
from lipyd.reader.peaks import peaks_sample_id_method

class TestSampleId(object):
    
    def test_plate_sample_id_processor(self):
        
        sip = sampleattrs.plate_sample_id_processor()
        
        a11_from_string = sip('A11')
        a11_from_tuple  = sip(('A', 11))
        
        assert a11_from_string == ('A', 11)
        assert a11_from_tuple  == ('A', 11)
    
    def test_generic_sample_id_processor(self):
        
        def _method(something):
            
            if isinstance(something, basestring):
                
                retime = re.compile(r'([0-9]+)\s?(s|min|h)')
                m = retime.search(something)
                
                if m:
                    
                    m = m.groups()
                    
                    return (int(m[0]), m[1])
            
            return something, None

        sip = sampleattrs.sample_id_processor(_method, 'time', 'unit')

        min10 = sip('10 min')
        h6    = sip('6h')
        other = sip('foobar')
        
        assert min10 == (10, 'min')
        assert h6    == (6,  'h')
        assert other == ('foobar', None)
    
    def test_passthrough_sample_id_processor(self):
        
        sip = sampleattrs.sample_id_processor()
        
        foobar = sip('foobar')
        
        assert hasattr(foobar, 'sample_id')
        assert foobar.sample_id == 'foobar'


class TestSampleAttrs(object):
    
    def test_sample_attrs_default(self):
        
        np.random.seed(123)
        sa = sampleattrs.SampleAttrs()
        
        assert sa.sample_id.sample_id == 'nccgrtkzwb'
    
    def test_sample_attrs_from_attrs(self):
        
        attrs = {
            'label': {
                'sample_id': ('A', 11),
            }
        }
        
        sa = sampleattrs.SampleAttrs(
            sample_id = peaks_sample_id_method,
            attrs = attrs,
            proc = sampleattrs.plate_sample_id_processor(),
        )
        
        assert sa.sample_id == ('A', 11)
    
    def test_sample_attrs_from_string(self):
        
        sa = sampleattrs.SampleAttrs(
            sample_id = ('A', 11),
            proc = sampleattrs.plate_sample_id_processor(),
        )
        
        assert sa.sample_id == ('A', 11)
    
    def test_sampleset_attrs_default(self):
        
        np.random.seed(123)
        ssa = sampleattrs.SampleSetAttrs(length = 9)
        
        assert len(ssa.sample_index_to_id) == 9
        assert ssa.sample_id_to_index[('hlwhblfxzs',)] == 7
    
    def test_sampleset_from_sample_ids_proc(self):
        
        ssa = sampleattrs.SampleSetAttrs(
            sample_ids = ['A9', 'A10', 'A11', 'A12', 'B1'],
            proc = sampleattrs.plate_sample_id_processor()
        )
        
        assert ssa.sample_index_to_id[-1].row == 'B'
        assert ssa.sample_id_to_index[('A', 12)] == 3
        assert len(ssa.sample_index_to_id) == 5
    
    def test_sampleset_from_attrs(self):
        
        ssa = sampleattrs.SampleSetAttrs(
            sample_ids = peaks_sample_id_method,
            attrs = (
                {'label': {'sample_id': 'G2',}},
                {'label': {'sample_id': 'G3',}},
                {'label': {'sample_id': 'G4',}},
            ),
            proc = sampleattrs.plate_sample_id_processor()
        )
        
        assert ssa.attrs[0].attrs['label']['sample_id'] == 'G2'
        assert ssa.sample_index_to_id[-1].col == 4
        assert len(ssa.sample_id_to_index) == 3
        assert ('G', 2) in ssa.sample_id_to_index
    
    def test_sampleset_attrs_argsort(self):
        
        ssa = sampleattrs.SampleSetAttrs(
            sample_ids = ['A9', 'A10', 'A11', 'A12', 'B1'],
            proc = sampleattrs.plate_sample_id_processor()
        )
        
        idx = ssa.argsort_by_sample_id(['A10', 'A11', 'A12', 'B1', 'A9'])
        
        assert np.all(idx == np.array([1, 2, 3, 4, 0]))
    
    def test_sampleset_attrs_sort(self):
        
        ssa = sampleattrs.SampleSetAttrs(
            sample_ids = ['A9', 'A10', 'A11', 'A12', 'B1'],
            proc = sampleattrs.plate_sample_id_processor(),
            attrs = (
                {'fr': 'a9'},
                {'fr': 'a10'},
                {'fr': 'a11'},
                {'fr': 'a12'},
                {'fr': 'b1'},
            ),
        )
        
        idx = ssa.argsort_by_sample_id(['A10', 'A11', 'A12', 'B1', 'A9'])
        
        ssa.sort_by_index(idx)
        
        assert ssa.sample_id_to_index[('A', 12)] == 2
        assert ssa.sample_index_to_id[0] == ('A', 10)
        assert ssa.attrs[-1].attrs['fr'] == 'a9'
    
    def test_sampleset_attrs_sort_by_id(self):
        
        ssa = sampleattrs.SampleSetAttrs(
            sample_ids = ['A9', 'A10', 'A11', 'A12', 'B1'],
            proc = sampleattrs.plate_sample_id_processor(),
        )
        
        ssa.sort_by_sample_id(['A10', 'A11', 'A12', 'B1', 'A9'])
        
        assert ssa.sample_id_to_index[('A', 12)] == 2
        assert ssa.sample_index_to_id[0] == ('A', 10)
    
    def test_sampleset_attrs_sort_by_other(self):
        
        ssa0 = sampleattrs.SampleSetAttrs(
            sample_ids = ['A9', 'A10', 'A11', 'A12', 'B1'],
            proc = sampleattrs.plate_sample_id_processor(),
        )
        
        ssa0.sort_by_sample_id(['A10', 'A11', 'A12', 'B1', 'A9'])
        
        ssa1 = sampleattrs.SampleSetAttrs(
            sample_ids = ['A9', 'A10', 'A11', 'A12', 'B1'],
            proc = sampleattrs.plate_sample_id_processor(),
        )
        
        ssa1.sort_by(ssa0)
        
        assert np.all(ssa0.sample_index_to_id == ssa1.sample_index_to_id)
        assert ssa1.sample_id_to_index[('A', 12)] == 2
        assert ssa1.sample_index_to_id[0] == ('A', 10)
        assert all(
            ssa0.sample_id_to_index[s] == ssa1.sample_id_to_index[s]
            for s in ssa0.sample_index_to_id
        )
    
    def test_sec_profile_1(self):
        
        peakspath = settings.get('peaks_gltpd1_invitro')
        secpath = settings.get('sec_gltpd1_invitro')
        
        reader = sample.SampleReader(
            input_type = 'peaks',
            fname = peakspath,
        )
        
        samples = reader.get_sampleset(
            sampleset_args = {
                'sample_id_proc': sampleattrs.plate_sample_id_processor(),
            }
        )
        
        secprofile = sampleattrs.SECProfile(
            sec_path = secpath,
            samples = samples,
        )
        
        assert secprofile.numof_samples == samples.numof_samples
        assert np.all(
            secprofile.attrs.sample_index_to_id ==
            samples.attrs.sample_index_to_id
        )
        assert id(samples) in secprofile._sample_data
        assert id(secprofile) in samples._sample_data
        assert secprofile.profile.max() - 143.40397368421048 < 0.0001
    
    def test_sec_profile_2(self):
        
        peakspath = settings.get('peaks_gltpd1_invivo')
        secpath = settings.get('sec_gltpd1_invivo')
        
        reader = sample.SampleReader(
            input_type = 'peaks',
            fname = peakspath,
        )
        
        samples = reader.get_sampleset(
            sampleset_args = {
                'sample_id_proc': sampleattrs.plate_sample_id_processor(),
            }
        )
        
        secprofile = sampleattrs.SECProfile(
            sec_path = secpath,
            samples = samples,
            start_volume = 1.2,
            offsets = (0.015, 0.045),
            start_col = 9,
            start_row = 'A',
            length = samples.numof_samples,
        )
        
        assert secprofile.numof_samples == samples.numof_samples
        assert np.all(
            secprofile.attrs.sample_index_to_id ==
            samples.attrs.sample_index_to_id
        )
        assert id(samples) in secprofile._sample_data
        assert id(secprofile) in samples._sample_data
        assert secprofile.profile015.max() - 20.92373913043478 < 0.0001
        assert secprofile.profile045.argmax() == 2
        assert secprofile.profile015.argmax() == 3
        assert secprofile.profiles == {'profile045', 'profile015'}
    
    def test_protein_containing_samples(self):
        
        peakspath = settings.get('peaks_gltpd1_invivo')
        secpath = settings.get('sec_gltpd1_invivo')
        
        reader = sample.SampleReader(
            input_type = 'peaks',
            fname = peakspath,
        )
        
        samples = reader.get_sampleset(
            sampleset_args = {
                'sample_id_proc': sampleattrs.plate_sample_id_processor(),
            }
        )
        
        secprofile = sampleattrs.SECProfile(
            sec_path = secpath,
            samples = samples,
            start_volume = 1.2,
            offsets = (0.015, 0.045),
            start_col = 9,
            start_row = 'A',
            length = samples.numof_samples,
        )
        
        pcs = secprofile.protein_containing_samples()
        
        assert np.all(
            pcs.selection ==
            np.array([False, False,  True,  True, False])
        )
        
        pcs = secprofile.protein_containing_samples(manual = ['A10', 'A11'])
        
        assert np.all(
            pcs.selection ==
            np.array([False,  True,  True, False, False])
        )
        
        pcs = secprofile.protein_containing_samples(exclude = ['A12'])
        
        assert np.all(
            pcs.selection ==
            np.array([False, False,  True, False, False])
        )
    
    def test_selection_logical_not(self):
        
        sel = sampleattrs.SampleSelection(
            selection = [False, True, True, False, False],
            sample_ids = ['A9', 'A10', 'A11', 'A12', 'B1'],
            sample_id_proc = sampleattrs.plate_sample_id_processor()
        )
        
        nsel = sel.logical_not()
        
        assert np.all(
            nsel.selection ==
            np.array([ True, False, False,  True,  True])
        )
