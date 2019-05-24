#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `lipyd` python module
#
#  Copyright (c) 2015-2019 - EMBL
#
#  File author(s):
#  Dénes Türei (turei.denes@gmail.com)
#  Igor Bulanov
#
#  Distributed under the GNU GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: https://saezlab.github.io/lipyd
#

import imp

import lipyd.session as session
import lipyd.msproc as msproc
import lipyd.sample as sample


class Main(session.Logger):
    
    
    def __init__(
        self,
        input_path,
        ionmode = None,
        sample_id_method = None,
        preprocess_args = None,
        ms2_param = None,
        **kwargs
    ):
        
        if not hasattr(self, '_log_name'):
            
            session.Logger.__init__(self, name = 'main')
        
        self.input_path = input_path
        self.ionmode = ionmode
        self.sample_id_method = None
        self.preprocess_args = preprocess_args or {}
        self.ms2_param = ms2_param
    
    
    def reload(self):
        
        modname = self.__class__.__module__
        mod = __import__(modname, fromlist = [modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)
    
    
    def main(self):
        
        self.preprocess()
        self.identification()
    
    
    def preprocess(self):
        
        self._log('Starting LC MS/MS data preprocessing.')
        
        self.preproc = msproc.MSPreprocess(
            input_path = self.input_path,
            ionmode = self.ionmode,
            sample_id_method = self.sample_id_method,
            **self.preprocess_args
        )
        self.preproc.main()
        
        self._log('Preprocessing finished.')
        
        self.init_sampleset()
    
    
    def init_sampleset(self):
        
        self._log('Creating SampleSet object.')
        
        sampleset_args = self.preproc.extractor.get_sampleset()
        feature_attrs = self.preproc.extractor.get_attributes()
        feattrs = sample.FeatureAttributes(**feature_attrs)
        sampleset_args['feature_attrs'] = feattrs
        sampleset_args['sorter'] = feattrs.sorter
        sampleset_args['ms2_param'] = self.ms2_param
        sampleset_args['ionmode'] = self.ionmode
        self.samples = sample.SampleSet(**sampleset_args)
        
        self._log(
            'SampleSet object created with %u features and %u samples.' % (
                len(self.samples),
                self.samples.numof_samples,
            )
        )
    
    
    def identification(self):
        
        self.samples.database_lookup()
        self.samples.set_ms2_param()
        self.samples.set_ms2_sources()
        self.samples.ms2_analysis()
