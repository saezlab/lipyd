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
#  Website: http://denes.omnipathdb.org/

import imp
import pandas as pd
import numpy as np
import lipyd.recalibration as recalibration
import lipyd.experiment as experiment
import lipyd.lookup as lookup
import csv
import collections


StandardCompoundsBase = collections.namedtuple(
    'StandardCompoundsBase',
    ['name', 'mass', 'formula', 'ionmode', 'adducts'] 
    )

filename = '/home/igor/Documents/Data/Alb_Lar/csv_neg/CERT_neg.csv'
mgfdir = '/home/igor/Documents/Data/Alb_Lar/mgf/190715_Olive_LvEk_AFO_CERT_infected_neg.mgf'


class StandardCompounds(object):

    def __init__(
        self,
        input_path = None,
        standards = None,
        ionmode = None,
        output_filename = None,
        mgfdir = None,
        input_ext = 'mzML'
    ):

        self.output_filename = output_filename
        self.input_path = input_path
        self.standards = standards or [] # [786.113, 810.025]
        self.ionmode = ionmode
        self.mgfdir = mgfdir
        self.input_ext = input_ext

    def mgf_match_method(path, attrs):
        
        return attrs.sample_id.sample_id in path

    def preproc(self):

        self.m = experiment.Experiment(
            nput_path = self.input_path,
            preprocess_args = {
                'input_filter': lambda n: self.ionmode in n,
                'input_ext': self.input_ext,
            },
            ionmode = self.ionmode,
            ms2_param = {
                'mgfdir': self.mgfdir,
                'mgf_match_method': self.mgf_match_method,
                'check_rt': True,
            }
        )
        
        self.m.preprocess()

    def lookup_standards(self):

        mztheo = []

        dfo = self.m.preproc.extractor.dataframe
           
        dfo.sort_values(by = 'mz', inplace = True)
        dfo.reset_index(drop = True, inplace = True)

        mzo = np.array(dfo['mz'])  

        for i in self.standards:
            a = lookup.find(mzo, i, t = 50)
            mztheo.append(a)

    def calc_ppm(self):

        ppms = ((mzo / mztheo) - 1) * 10e6  

    def run(self):

        self.preproc()
        self.lookup_standards()
        self.calc_ppm()

    def reload(self):
        """ """
        
        modname = self.__class__.__module__
        mod = __import__(modname, fromlist = [modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)