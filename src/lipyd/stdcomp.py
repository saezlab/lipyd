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


import lipyd.recalibration as recalibration
import lipyd.experiment as experiment
import lipyd.lookup as lookup
import csv
import collections


StandardCompoundsBase = collections.namedtuple(
	'StandardCompoundsBase',
	['name', 'mass', 'formula', 'ionmode', 'adducts'] 
	)

class StandardCompounds(collections.namedtuple(
	'StandardCompoundsBase',
	[
		'name', 
		'mass', 
		'formula', 
		'ionmode', 
		'adducts'
	] 
	)):

    def __init__(self,
            input_path = None,
            standards = None,
            ionmode = None,
            output_filename = None
        ):

        self.output_filename = output_filename
        self.input_path = input_path
        self.standards = standards or {}
        self.ionmode = ionmode

    def preprocess(self, ionmode):

        m = experiment.Experiment(
        input_path = (
            self.input_path
        ),
        preprocess_args = {
            'input_filter': lambda n: self.ionmode in n,
            'input_ext': 'mzML',
        },
        ionmode = self.ionmode,
        )

        m.preprocess()

        m.export(self.output_filename)

    def lookup_standards(self):

        pairs = set()

        df = csv.read(self.output_filename)

        for k,v in self.standards:
            a = lookup.find(k, mz, t = 50)
            pairs.add((k,a))

    

